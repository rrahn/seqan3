// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <vector>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_element_affine_cell.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

// ----------------------------------------------------------------------------
// global types
// ----------------------------------------------------------------------------

using namespace seqan3::detail;
using namespace seqan3;

using score_t = int16_t;
using simd_score_t = simd_type_t<score_t>;
using alphabet_t = dna4;
using score_pair_t = std::pair<simd_score_t, simd_score_t>;
using cell_score_t = score_pair_t;

// ----------------------------------------------------------------------------
// sparse matrix
// ----------------------------------------------------------------------------

// Can we optimise the matrix?
template <typename value_t>
class sparse_matrix
{
private:

    // using scalar_t = typename simd_traits<value_t>::scalar_type;
    // using array_t = std::array<scalar_t, simd_traits<value_t>::length>;
    using element_t = std::pair<value_t, value_t>;
    using column_t = std::vector<element_t, aligned_allocator<element_t, simd_traits<value_t>::max_length>>;

    class iterator;

    column_t matrix_data;
    size_t col_dim{};
public:

    template <typename idx_t>
    void resize(column_index_type<idx_t> col_dim, row_index_type<idx_t> row_dim)
    {
        matrix_data.resize(row_dim.get());
        this->col_dim = col_dim.get();
    }

    auto begin()
    {
        return iterator(matrix_data, 0u);
    }

    auto end()
    {
        return iterator(matrix_data, col_dim);
    }
};

template <typename value_t>
class sparse_matrix<value_t>::iterator
{
private:

    std::span<element_t> column{};
    size_t current_column_id{};
public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = std::span<element_t>;
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    iterator() noexcept = default; //!< Defaulted.

    explicit iterator(column_t & matrix_data, size_t const initial_column_id) noexcept :
        column{matrix_data},
        current_column_id{initial_column_id}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        return column;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    iterator & operator++()
    {
        ++current_column_id;
        return *this;
    }

    //!\brief Move `this` to the next column.
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Tests whether `lhs == rhs`.
    bool operator==(iterator const & rhs) const noexcept
    {
        return current_column_id == rhs.current_column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    bool operator!=(iterator const & rhs) const noexcept
    {
        return !(*this == rhs);
    }
    //!\}
};

// ----------------------------------------------------------------------------
// compute matrix
// ----------------------------------------------------------------------------

template <typename cell_t>
inline auto compute_cell(simd_score_t diagonal_score,
                         cell_t cell,
                         simd_score_t score,
                         simd_score_t & vertical_score,
                         simd_score_t gap_extension_score,
                         simd_score_t gap_open_score)
{
    // using byte_t = typename cell_t::first_type;

    diagonal_score += score;
    // simd_score_t horizontal_score = cell.second;
    // simd_score_t vertical_score = cell.vertical_score();
    // simd_score_t horizontal_score = reinterpret_cast<simd_score_t &>(cell.second);

    diagonal_score = (diagonal_score < vertical_score) ? vertical_score : diagonal_score;
    diagonal_score = (diagonal_score < cell.second) ? cell.second : diagonal_score;

    score = diagonal_score + gap_open_score;
    vertical_score += gap_extension_score;
    cell.second += gap_extension_score;

    // store the vertical_score and cell.second value in the next path
    vertical_score = (vertical_score < score) ? score : vertical_score;
    cell.second = (cell.second < score) ? score : cell.second;

    return cell_score_t{diagonal_score, cell.second};
    // return cell_t{reinterpret_cast<byte_t &>(diagonal_score), reinterpret_cast<byte_t &>(cell.second)};
}

template <typename seq1_t, typename seq2_t, typename matrix_t, typename scoring_scheme_t >
inline auto compute_matrix(seq1_t && seq1,
                           seq2_t && seq2,
                           matrix_t & matrix,
                           scoring_scheme_t const & scheme,
                           score_t gap_extension,
                           score_t gap_open)
{
    simd_score_t gap_extension_score = simd::fill<simd_score_t>(gap_extension);
    simd_score_t gap_open_score = simd::fill<simd_score_t>(gap_open);

    auto matrix_it = matrix.begin(); // Get the first iterator.
    auto && column = *matrix_it;

    using pure_data_t = std::ranges::range_value_t<decltype(column)>;

    // ------------------------------------------------------------------------
    // initilaise the first column

    simd_score_t vertical_score{};
    for (auto col_it = column.begin(); col_it != column.end(); ++col_it)
        *col_it = pure_data_t{}; // everything is 0

    // ------------------------------------------------------------------------
    // compute remaining matrix

    for (auto seq1_it = (++matrix_it, seq1.begin()); matrix_it != matrix.end(); ++matrix_it, ++seq1_it)
    {
        column = *matrix_it;
        auto col_it = column.begin(); // iterator to the cell in the first row of the current column
        // initialise first value
        simd_score_t diagonal = ((*col_it).first); // We will use at the beginning our current score matrix implementation
        *col_it = pure_data_t{}; // we use zero intialisation again.
        vertical_score = simd_score_t{}; // initialise vertical score.

        // Move down remaining column.
        for (auto seq2_it = (++col_it, seq2.begin()); col_it != column.end(); ++col_it, ++seq2_it)
        {
            simd_score_t next_diagonal = ((*col_it).first);
            simd_score_t score = scheme.score(*seq1_it, *seq2_it);
            *col_it = compute_cell(diagonal, *col_it, score, vertical_score, gap_extension_score, gap_open_score);
            diagonal = next_diagonal;
        }
    }

    auto it = std::next(column.begin(), std::ranges::distance(column) - 1);
    return ((*it).first);
}

// ----------------------------------------------------------------------------
// transform sequence batch to simd.
// ----------------------------------------------------------------------------

template <typename seq_t, typename padding_symbol_t>
inline auto to_simd(seq_t && seq_batch, padding_symbol_t const & padding_symbol)
{
    using simd_seq_t = std::vector<simd_score_t, aligned_allocator<simd_score_t, alignof(simd_score_t)>>;
    simd_seq_t simd_sequence{};

    for (auto && simd_vector_chunk : seq_batch | views::to_simd<simd_score_t>(padding_symbol))
        std::ranges::move(simd_vector_chunk, std::cpp20::back_inserter(simd_sequence));

    return simd_sequence;
}

void simd_score_bench(benchmark::State & state)
{
    // ----------------------------------------------------------------------------
    // Create scoring scheme
    // ----------------------------------------------------------------------------

    using scoring_scheme_t = simd_match_mismatch_scoring_scheme<simd_score_t, alphabet_t, align_cfg::method_global>;

    scoring_scheme_t scoring_scheme{nucleotide_scoring_scheme{match_score<score_t>{4}, mismatch_score<score_t>{-5}}};

    // ----------------------------------------------------------------------------
    // Generate and transform sequences
    // ----------------------------------------------------------------------------

    // One batch and how often can be compute the alignment.
    auto data = seqan3::test::generate_sequence_pairs<seqan3::dna4>(500, simd_traits<simd_score_t>::length);

    auto simd_seq1 = to_simd(data | views::get<0>, scoring_scheme_t::padding_symbol);
    auto simd_seq2 = to_simd(data | views::get<1>, scoring_scheme_t::padding_symbol);

    // ----------------------------------------------------------------------------
    // Enable score matrix
    // ----------------------------------------------------------------------------

    // using matrix_element_t = alignment_matrix_element_affine_cell<simd_score_t>;
    // using matrix_t = score_matrix_single_column<matrix_element_t>;
    using matrix_t = sparse_matrix<simd_score_t>;

    matrix_t matrix{};
    matrix.resize(column_index_type{simd_seq1.size()}, row_index_type{simd_seq2.size()});

    // ----------------------------------------------------------------------------
    // Compute alignment
    // ----------------------------------------------------------------------------

    simd_score_t final_score{};
    for (auto _ : state)
        final_score = compute_matrix(simd_seq1, simd_seq2, matrix, scoring_scheme, -1, -11);

    benchmark::DoNotOptimize(final_score);

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(data, configuration{align_cfg::method_global{}});
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(simd_score_bench);

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
