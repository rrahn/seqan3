// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/type_inspection.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/zip.hpp>

#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The configuration type; must be of type seqan3::configuration.
 * \tparam policies_t Variadic template argument for the different policies of this alignment algorithm.
 *
 * \details
 *
 * ### Configuration
 *
 * The first template argument is the type of the alignment configuration. The alignment configuration was used to
 * configure the `alignment algorithm type` within the seqan3::detail::alignment_configurator.
 * The algorithm computes a column based dynamic programming matrix given two sequences.
 * After the computation a user defined callback function is invoked with the computed seqan3::alignment_result.
 */
template <typename alignment_configuration_t, typename ...policies_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class pairwise_alignment_algorithm : protected policies_t...
{
protected:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using traits_type = configuration_traits_t<alignment_configuration_t>;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;
    //!\brief The configured alignment result type.
    using alignment_result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");


public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm() = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the base policies of the alignment algorithm.
     */
    pairwise_alignment_algorithm(alignment_configuration_t const & config) : policies_t(config)...
    {}
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given range over indexed sequence pairs.
     * \tparam indexed_sequence_pairs_t The type of indexed_sequence_pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \tparam callback_t The type of the callback function that is called with the alignment result; must model
    *                    std::invocable with seqan3::alignment_result as argument.
     *
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     * \param[in] callback The callback function to be invoked with each computed alignment result.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment for each
     * sequence pair. The space and runtime complexities depend on the selected configurations (see below).
     * For every computed alignment the given callback is invoked with the respective alignment result.
     *
     * ### Exception
     *
     * Strong exception guarantee. Might throw std::bad_alloc or seqan3::invalid_alignment_configuration.
     *
     * ### Thread-safety
     *
     * Calls to this functions in a concurrent environment are not thread safe. Instead use a copy of the alignment
     * algorithm type.
     *
     * ### Complexity
     *
     * The following table lists the runtime and space complexities for the banded and unbanded algorithm dependent
     * on the configured seqan3::align_cfg::result per sequence pair.
     * Let `n` be the length of the first sequence, `m` be the length of the second sequence and `k` be the size of
     * the band.
     *
     * |                        | unbanded         | banded            |
     * |:----------------------:|:----------------:|:-----------------:|
     * |runtime                 |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (score only)      |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (end positions)   |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (begin positions) |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (alignment)       |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     */
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires std::invocable<callback_t, alignment_result_type>
    //!\endcond
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_type>::type;
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));

            result_value_t res{};
            res.id = idx;
            res.score = this->optimal_score;
            callback(alignment_result_type{res});
        }
    }

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires traits_type::is_vectorised && std::invocable<callback_t, alignment_result_type>
    //!\endcond
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_type>::type;
        using simd_collection_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;

        // Extract the batch of sequences for the first and the second sequence.
        auto seq1_collection = indexed_sequence_pairs | views::get<0> | views::get<0>;
        auto seq2_collection = indexed_sequence_pairs | views::get<0> | views::get<1>;

        this->initialise_tracker(seq1_collection, seq2_collection);

        // Convert batch of sequences to sequence of simd vectors.
        thread_local simd_collection_t simd_seq1_collection{};
        thread_local simd_collection_t simd_seq2_collection{};

        convert_batch_of_sequences_to_simd_vector(simd_seq1_collection, seq1_collection, traits_type::padding_symbol);
        convert_batch_of_sequences_to_simd_vector(simd_seq2_collection, seq2_collection, traits_type::padding_symbol);

        compute_matrix(simd_seq1_collection, simd_seq2_collection);

        size_t index = 0;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            (void) sequence_pair;
            result_value_t res{};
            res.id = idx;

            if constexpr (traits_type::result_type_rank >= 0)
            {
                res.score = this->optimal_score[index] -
                            (this->padding_offsets[index] * this->m_scoring_scheme.max_match_score());
            }

            callback(alignment_result_type{res});
            ++index;
        }
    }

protected:
    /*!\brief Converts a batch of sequences to a sequence of simd vectors.
     * \tparam simd_sequence_t The type of the simd sequence; must model std::ranges::output_range for the `score_type`.
     * \tparam sequence_collection_t The type of the collection containing the sequences; must model
     *                               std::ranges::forward_range.
     * \tparam padding_symbol_t The type of the padding symbol.
     *
     * \param[out] simd_sequence The transformed simd sequence.
     * \param[in] sequences The batch of sequences to transform.
     * \param[in] padding_symbol The symbol that should be appended during the transformation/
     *
     * \details
     *
     * Expects that the size of the collection is less or equal than the number of alignments that can be computed
     * within one simd vector (simd_traits\<score_type\>\::length).
     * Applies an Array-of-Structures (AoS) to Structure-of-Arrays (SoA) transformation by storing one
     * column of the collection as a simd vector. The resulting simd sequence has the size of the longest sequence in
     * the collection. For all sequences with a smaller size the padding symbol will be appended during the simd
     * transformation to fill up the remaining size difference.
     */
    template <typename simd_sequence_t,
              std::ranges::forward_range sequence_collection_t,
              arithmetic padding_symbol_t>
    //!\cond
        requires std::ranges::output_range<simd_sequence_t, score_type>
    //!\endcond
    void convert_batch_of_sequences_to_simd_vector(simd_sequence_t & simd_sequence,
                                                   sequence_collection_t & sequences,
                                                   padding_symbol_t const & padding_symbol)
    {
        assert(static_cast<size_t>(std::ranges::distance(sequences)) <= traits_type::alignments_per_vector);

        simd_sequence.clear();
        for (auto && simd_vector_chunk : sequences | views::to_simd<score_type>(padding_symbol))
            std::ranges::move(simd_vector_chunk, std::ranges::back_inserter(simd_sequence));
    }

    /*!\brief Compute the actual alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void compute_matrix(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------
        using matrix_index_t = typename traits_type::matrix_index_type;

        this->reset_optimum(); // Reset the tracker for the new alignment computation.

        thread_local score_matrix_single_column<score_type> local_score_matrix{};
        coordinate_matrix<matrix_index_t> local_index_matrix{};

        size_t number_of_columns = std::ranges::distance(sequence1) + 1;
        size_t number_of_rows = std::ranges::distance(sequence2) + 1;

        local_score_matrix.resize(column_index_type{number_of_columns}, row_index_type{number_of_rows});
        local_index_matrix.resize(column_index_type{number_of_columns}, row_index_type{number_of_rows});

        auto alignment_matrix_it = local_score_matrix.begin();
        auto indexed_matrix_it = local_index_matrix.begin();

        initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        size_t seq1_size = sequence1.size() - 1;
        auto it_seq1 = sequence1.begin();
        for (size_t idx = 0; idx < seq1_size; ++idx, ++it_seq1)
            compute_column(*++alignment_matrix_it, *++indexed_matrix_it, *it_seq1, sequence2);

        compute_last_column(*++alignment_matrix_it, *++indexed_matrix_it, *it_seq1, sequence2);

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        // auto && alignment_column = *alignment_matrix_it;
        // auto && cell_index_column = *indexed_matrix_it;

        // auto alignment_column_it = alignment_column.begin();
        // auto cell_index_column_it = cell_index_column.begin();

        // this->track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        // for ([[maybe_unused]] auto && unused : sequence2)
        //     this->track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        // this->track_final_cell(*alignment_column_it, *cell_index_column_it);
    }

    /*!\brief Initialise the first column of the alignment matrix.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::input_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::input_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::input_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] sequence2 The second sequence used to determine the size of the column.
     *
     * \details
     *
     * The first column of the alignment matrix does not require any character comparisons of the sequences that
     * shall be aligned. The second sequence is thus only needed to determine the size of the column.
     * The computation of the column is split into three phases: the initialisation phase, the iteration phase, and
     * the final phase. In the initialisation phase the first cell of the column is computed and in the iteration
     * phase all remaining cells are computed. In the final phase the last cell is possibly evaluated for a new
     * alignment optimum.
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              std::ranges::input_range sequence2_t>
    void initialise_column(alignment_column_t && alignment_column,
                           cell_index_column_t && cell_index_column,
                           sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto first_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();
        *first_column_it = this->track_cell(this->initialise_origin_cell(), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for ([[maybe_unused]] auto const & unused : sequence2)
        {
            ++first_column_it;
            *first_column_it = this->track_cell(this->initialise_first_column_cell(*first_column_it),
                                                *++cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell of initial column
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*first_column_it, *cell_index_column_it);

        // auto alignment_column_it = alignment_column.begin();
        // std::cout << "<D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";

        // for ([[maybe_unused]] auto const & unused : sequence2)
        // {
        //     ++alignment_column_it;
        //     std::cout << " <D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";
        // }
        // std::cout << "\n";
    }

    /*!\brief Initialise any column of the alignment matrix except the first one.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::input_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::input_range.
     * \tparam sequence1_value_t The value type of sequence1; must model seqan3::semialphabet.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::input_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] sequence1_value The current symbol of sequence1.
     * \param[in] sequence2 The second sequence to align against `sequence1_value`.
     *
     * \details
     *
     * Computes the alignment for the given alignment matrix column. The function splits the computation of the column
     * into three phases: the initialisation phase, the iteration phase, and the final phase. In the initialisation
     * phase the first cell of the column is computed and in the iteration phase all remaining cells are computed.
     * In the final phase the last cell is possibly evaluated for a new alignment optimum.
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              typename sequence1_value_t,
              std::ranges::input_range sequence2_t>
    //!\cond
        requires semialphabet<sequence1_value_t> || simd_concept<sequence1_value_t>
    //!\endcond
    void compute_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        sequence1_value_t const & sequence1_value,
                        sequence2_t && sequence2)
    {
        using score_type = typename traits_type::score_type;

        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        score_type diagonal = cell.optimal_score();
        *alignment_column_it = this->track_cell(this->initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto const & sequence2_value : sequence2)
        {
            auto cell = *++alignment_column_it;
            score_type next_diagonal = cell.optimal_score();
            *alignment_column_it = this->track_cell(
                this->compute_inner_cell(diagonal, cell, this->m_scoring_scheme.score(sequence1_value, sequence2_value)),
                *++cell_index_column_it);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*alignment_column_it, *cell_index_column_it);

        // alignment_column_it = alignment_column.begin();
        // std::cout << "<D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";

        // for ([[maybe_unused]] auto const & unused : sequence2)
        // {
        //     ++alignment_column_it;
        //     std::cout << " <D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";
        // }
        // std::cout << "\n";
    }

    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              typename sequence1_value_t,
              std::ranges::input_range sequence2_t>
    //!\cond
        requires semialphabet<sequence1_value_t> || simd_concept<sequence1_value_t>
    //!\endcond
    void compute_last_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        sequence1_value_t const & sequence1_value,
                        sequence2_t && sequence2)
    {
        using score_type = typename traits_type::score_type;

        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        score_type diagonal = cell.optimal_score();
        *alignment_column_it = this->track_last_column_cell(this->initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto const & sequence2_value : sequence2)
        {
            auto cell = *++alignment_column_it;
            score_type next_diagonal = cell.optimal_score();
            *alignment_column_it = this->track_last_column_cell(
                this->compute_inner_cell(diagonal, cell, this->m_scoring_scheme.score(sequence1_value, sequence2_value)),
                *++cell_index_column_it);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*alignment_column_it, *cell_index_column_it);
        this->track_final_cell(*alignment_column_it, *cell_index_column_it);

        // alignment_column_it = alignment_column.begin();
        // std::cout << "<D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";

        // for ([[maybe_unused]] auto const & unused : sequence2)
        // {
        //     ++alignment_column_it;
        //     std::cout << " <D:" << (int)(*alignment_column_it).optimal_score()[0] << " H:" << (int)(*alignment_column_it).horizontal_score()[0] << " V:" <<  (int)(*alignment_column_it).vertical_score()[0] << ">";
        // }
        // std::cout << "\n";
    }
};

} // namespace seqan3::detail
