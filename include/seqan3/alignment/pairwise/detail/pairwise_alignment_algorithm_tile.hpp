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
#include <chrono>

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
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/zip.hpp>

#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>

#include <seqan3/core/simd/debug_stream_simd.hpp>
#include <seqan3/core/debug_stream.hpp>

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
class pairwise_alignment_algorithm_tile :
    protected pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>
{
protected:
    //!\brief The type of the algorithm to inherit from.
    using base_algorithm_t = pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>;
    //Import base types.
    using typename base_algorithm_t::traits_type;
    using typename base_algorithm_t::score_type;
    using typename base_algorithm_t::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");

    using duration_type = typename std::chrono::high_resolution_clock::duration;

    // scoring_scheme_type scoring_scheme;
    // duration_type duration_prep{};
    // duration_type duration_algo{};
    // duration_type duration_res{};

    size_t tile_size{};
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm_tile() = default; //!< Defaulted.
    pairwise_alignment_algorithm_tile(pairwise_alignment_algorithm_tile const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm_tile(pairwise_alignment_algorithm_tile &&) = default; //!< Defaulted.
    pairwise_alignment_algorithm_tile & operator=(pairwise_alignment_algorithm_tile const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm_tile & operator=(pairwise_alignment_algorithm_tile &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm_tile()
    {
        // auto end_time = std::chrono::high_resolution_clock::now();
        // std::cout << "Preparation:    " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_prep).count() << " ms\n";
	    // std::cout << "Algorithm:      " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_algo).count() << " ms\n";
	    // std::cout << "Extract result: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_res).count()  << " ms\n";

    } //= default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the base policies of the alignment algorithm. Subsequently determines the tile size based on the
     * biggest match score, the lowest mismatch score and the selected recursion policy.
     */
    pairwise_alignment_algorithm_tile(alignment_configuration_t const & config) : base_algorithm_t(config)
    {
        tile_size = this->determine_tile_size(this->m_scoring_scheme.max_match_score(),
                                              this->m_scoring_scheme.min_mismatch_score());
    }
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
    using base_algorithm_t::operator();

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires traits_type::is_vectorised && std::invocable<callback_t, alignment_result_type>
    //!\endcond
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_type>::type;
        using simd_collection_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;
        using original_score_type = typename traits_type::original_score_type;

        // auto start = std::chrono::high_resolution_clock::now();
        // Extract the batch of sequences for the first and the second sequence.
        auto seq1_collection = indexed_sequence_pairs | views::get<0> | views::get<0>;
        auto seq2_collection = indexed_sequence_pairs | views::get<0> | views::get<1>;

        this->initialise_tracker(seq1_collection, seq2_collection, tile_size);

        // Convert batch of sequences to sequence of simd vectors.
        thread_local simd_collection_t simd_seq1_collection{};
        thread_local simd_collection_t simd_seq2_collection{};

        base_algorithm_t::convert_batch_of_sequences_to_simd_vector(simd_seq1_collection, seq1_collection,
                                                                    traits_type::padding_symbol);
        base_algorithm_t::convert_batch_of_sequences_to_simd_vector(simd_seq2_collection, seq2_collection,
                                                                    traits_type::padding_symbol);

        // duration_prep += std::chrono::high_resolution_clock::now() - start;
        // start = std::chrono::high_resolution_clock::now();

        compute_tile_matrix(simd_seq1_collection, simd_seq2_collection);

        // duration_algo += std::chrono::high_resolution_clock::now() - start;
        // start = std::chrono::high_resolution_clock::now();

        // std::cout << std::setw(20) << "score: "          << std::setw(6) << (int)this->optimal_score[0] << "\n";
        // std::cout << std::setw(20) << "row offset: "     << std::setw(6) << (int)this->score_row_offsets[0] << "\n";
        // std::cout << std::setw(20) << "col offset: "     << std::setw(6) << (int)this->score_column_offsets[0] << "\n";
        // std::cout << std::setw(20) << "padding offset: " << std::setw(6) << (int)this->padding_offsets[0] << "\n";

        size_t index = 0;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            (void) sequence_pair;
            result_value_t res{};
            res.id = idx;

            if constexpr (traits_type::result_type_rank >= 0)
            {
                // Based on the builder pattern?
                //TODO: This should be defined as a function over the policy.
                res.score = original_score_type{this->optimal_score[index]} +
                            this->score_row_offsets[index] +
                            this->score_column_offsets[index] -
                            (this->padding_offsets[index] * this->m_scoring_scheme.max_match_score());
            }

            callback(alignment_result_type{res});
            ++index;
        }
        // duration_res += std::chrono::high_resolution_clock::now() - start;
    }

protected:
    /*!\brief Compute the actual alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void compute_tile_matrix(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        // auto print_tile = [ ] (auto col, auto row)
        // {
        //     for (auto row_ : row)
        //         debug_stream << "<D:" << std::setw(4) << (int) row_.optimal_score()[0] << " "
        //                      <<  "H:" << std::setw(4) << (int) row_.horizontal_score()[0] << " "
        //                      <<  "V:" << std::setw(4) << (int) row_.vertical_score()[0] << "> ";

        //     debug_stream << "\n";
        //     for (auto col_ : col)
        //         debug_stream << "<D:" << std::setw(4) << (int) col_.optimal_score()[0] << " "
        //                      <<  "H:" << std::setw(4) << (int) col_.horizontal_score()[0] << " "
        //                      <<  "V:" << std::setw(4) << (int) col_.vertical_score()[0] << "> ";

        //     debug_stream << "\n\n";
        // };

        //----------------------------------------------------------------------
        // Initialise the tiles:
        //----------------------------------------------------------------------

        // We use the same index type as for the local tile matrices.
        using matrix_index_t = typename traits_type::matrix_index_type;

        // We need somehow the cell type.
        using cell_t = affine_cell_proxy<std::tuple<score_type, score_type, score_type>>;
        using cell_vector_t = std::vector<cell_t, aligned_allocator<cell_t, alignof(cell_t)>>;
        using offset_vector_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;

        thread_local cell_vector_t tile_row{};
        thread_local std::vector<cell_vector_t> tile_column_vector{};
        thread_local offset_vector_t tile_offsets{};

        size_t const tile_column_count = (sequence1.size() - 1 + tile_size) / tile_size;
        size_t const tile_row_count = (sequence2.size() - 1 + tile_size) / tile_size;

        // We won't use the first cell of the row.
        tile_row.resize(tile_size + 1);
        tile_column_vector.resize(tile_row_count, tile_row);
        tile_offsets.resize(tile_row_count);

        // We could use a coordinate matrix to reflect the tile ids.
        coordinate_matrix<matrix_index_t> tile_index_matrix{};
        tile_index_matrix.resize(column_index_type{tile_row_count}, row_index_type{tile_row_count});

        //----------------------------------------------------------------------
        // Initialise the first column:
        //----------------------------------------------------------------------

        size_t seq1_tile_position = 0;
        size_t seq2_tile_position = 0;

        auto sequence1_tile = sequence1 | views::slice(seq1_tile_position, seq1_tile_position + tile_size);
        auto tile_column_it = tile_column_vector.begin();
        auto tile_offset_it = tile_offsets.begin();
        auto index_matrix_it = tile_index_matrix.begin();

        auto index_column = *index_matrix_it;
        auto index_column_it = index_column.begin();
        this->initialise_first_tile(*tile_column_it, tile_row);

        // debug_stream << "\nBefore:\n";
        // print_tile(*tile_column_it, tile_row);

        this->set_tile_to_track(*tile_column_it, tile_row, *index_column_it);
        *tile_offset_it = score_type{};

        this->compute_matrix(sequence1_tile,
                             sequence2 | views::slice(seq2_tile_position, seq2_tile_position + tile_size));
        // this->track_tile();  // Special handling for searching the maximum.

        // debug_stream << "After:\n";
        // print_tile(*tile_column_it, tile_row);

        for (size_t tile_id = 1; tile_id < tile_row_count; ++tile_id)
        {
            seq2_tile_position += tile_size;
            ++index_column_it;
            ++tile_column_it;
            *++tile_offset_it = tile_row.front().optimal_score();
            this->initialise_first_column_tile(*tile_column_it, tile_row);

            // debug_stream << "Before:\n";
            // print_tile(*tile_column_it, tile_row);

            this->set_tile_to_track(*tile_column_it, tile_row, *index_column_it);
            this->compute_matrix(sequence1_tile,
                                 sequence2 | views::slice(seq2_tile_position, seq2_tile_position + tile_size));
            // this->track_inner_tile();
            // debug_stream << "After:\n";
            // print_tile(*tile_column_it, tile_row);
        }

        this->track_tile_last_row(tile_row,
                                  *index_column_it,
                                  this->compute_column_offset(*index_column_it, tile_size),
                                  tile_offsets);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        for (size_t tile_column_idx = 1; tile_column_idx < tile_column_count; ++tile_column_idx)
        {
            ++index_matrix_it;
            index_column = *index_matrix_it;
            index_column_it = index_column.begin();
            seq1_tile_position += tile_size;
            size_t seq2_tile_position = 0;

            auto sequence1_tile = sequence1 | views::slice(seq1_tile_position, seq1_tile_position + tile_size);
            tile_column_it = tile_column_vector.begin();
            tile_offset_it = tile_offsets.begin();

            // print_tile(*tile_column_it, tile_row);

            *tile_offset_it = score_type{}; // First offset is always 0, because we compute it later.
            this->initialise_first_row_tile(*tile_column_it, tile_row);

            // debug_stream << "Before:\n";
            // print_tile(*tile_column_it, tile_row);

            this->set_tile_to_track(*tile_column_it, tile_row, *index_column_it);
            this->compute_matrix(sequence1_tile,
                                 sequence2 | views::slice(seq2_tile_position, seq2_tile_position + tile_size));

            // debug_stream << "After:\n";
            // print_tile(*tile_column_it, tile_row);
            // this->track_inner_tile();
            ++tile_column_it;
            seq2_tile_position += tile_size;
            for (; tile_column_it != tile_column_vector.end(); ++tile_column_it)
            {
                ++index_column_it;
                *++tile_offset_it = tile_row.front().optimal_score();
                this->initialise_inner_tile(*tile_column_it, tile_row);

                // debug_stream << "Before:\n";
                // print_tile(*tile_column_it, tile_row);

                this->set_tile_to_track(*tile_column_it, tile_row, *index_column_it);
                this->compute_matrix(sequence1_tile,
                                     sequence2 | views::slice(seq2_tile_position, seq2_tile_position + tile_size));
                seq2_tile_position += tile_size;
                // debug_stream << "After:\n";
                // print_tile(*tile_column_it, tile_row);
                // this->track_inner_tile();
            }

            this->track_tile_last_row(tile_row,
                                      *index_column_it,
                                      this->compute_column_offset(*index_column_it, tile_size),
                                      tile_offsets);
        }

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        // We need to iterate over the last column except the last one, which was already covered before.
        index_column = *index_matrix_it;
        index_column_it = index_column.begin();
        for (size_t idx = 0; idx < tile_column_vector.size() - 1; ++idx, ++index_column_it)
            this->track_tile_last_column(tile_column_vector[idx],
                                         *index_column_it,
                                         this->compute_column_offset(*index_column_it, tile_size),
                                         tile_offsets);
    }
};

} // namespace seqan3::detail

//                              +14                     +4                      +4
// //      e,  A,  A,  C,  C,     e,  G,  G,  T,  T,      e,  A,  A,  C,  C,      e,  G,  G,  T,  T
// /*e*/  0 ,-11,-12,-13,-14,     0, -1, -2, -3, -4,      0, -1, -2, -3, -4,      0, -1, -2, -3, -4,
// /*A*/ -11,  4, -7, -8, -9,     5,  4,  3,  2,  1,      5,  4,  3,  2,  1,      5,  4,  3,  2,  1,
// /*C*/ -12, -7, -1, -3, -4,    10,  0, -1, -2, -3,      1,  0, -1,  7,  6,     10,  0, -1, -2, -3,
// /*G*/ -13, -8,-12, -6, -8,     6, 14,  4,  2,  1,      5,  4,  3,  2,  2,      6, 14,  4,  2,  1,
// /*T*/ -14, -9,-13,-15,-11,     3,  3,  9,  8,  6,     10,  0, -1, -2, -3,      1,  3,  9,  8,  6,

// /*e*/  0,  5,  1, -1,  3,      0,  0,  6,  5,  3,      0,-10,-11,-12,-13,      0,  2,  8,  7,  5,
// /*A*/ -1,  4,  9, -2, -3,     -6, -1, -5,  1,  0,     -3,  4, -6, -8, -9,      4,  5,  4,  5,  4,
// /*C*/ -2,  3, -1, 13,  2,     -1, -2, -3, -4, -4,     -7, -7, -1, -2, -4,      9,  2,  1,  0,  0,
// /*G*/ -3,  2, -2,  2,  8,      5,  3,  2, -8, -9,    -12, -8,-12, -6, -7,      6, 13,  6,  3,  2,
// /*T*/ -4,  1, -3,  1, -3,     -6,  0, -2,  6, -4,     -7, -9,-10,-11,-11,      2,  2, 10, 10,  7,

// /*e*/  0,  5,  1,  5,  1,      0,  6,  4, 12,  2,      0, -2, -3, -4, -4,      0,  2,  8,  8,  5,
// /*A*/ -1,  4,  9,  4,  0,     -1,  1,  1,  1,  7,      5,  4,  2, -8, -9,     -5, -1, -5,  1,  1,

// +0                                                                                                                   // +14                                                                                                                  // +4                                                                                                                  // +4
// <D:   0 H: -11 V: -11> <D: -11 H: -12 V: -22> <D: -12 H: -13 V: -23> <D: -13 H: -14 V: -24> <D: -14 H: -15 V: -25>   // <D: -14 H: -15 V: -25> <D: -15 H: -16 V: -26> <D: -16 H: -17 V: -27> <D: -17 H: -18 V: -28> <D:  -4 H:  -5 V: -15>   // <D:  -4 H:  -5 V: -15> <D:  -5 H:  -6 V: -16> <D:  -6 H:  -7 V: -17> <D:  -7 H:  -8 V: -18> <D:  -4 H:  -5 V: -15>  // <D:  -4 H:  -5 V: -15> <D:  -5 H:  -6 V: -16> <D:  -6 H:  -7 V: -17> <D:  -7 H:  -8 V: -18> <D:  -4 H:  -5 V: -15>
// <D: -11 H: -22 V: -12>                                                                      <D:  -9 H: -10 V: -20>   // <D:  -9 H: -10 V: -20>                                                                      <D:   1 H:   0 V: -10>   // <D:   1 H:   0 V: -10>                                                                      <D:   1 H:   0 V: -10>  // <D:   1 H:   0 V: -10>                                                                      <D:   1 H:   0 V: -10>
// <D: -12 H: -23 V: -13>                                                                      <D:  -4 H: -14 V: -15>   // <D:  -4 H: -14 V: -15>                                                                      <D:  -3 H:  -4 V: -11>   // <D:  -3 H:  -4 V: -11>                                                                      <D:   6 H:  -4 V:  -5>  // <D:   6 H:  -4 V:  -5>                                                                      <D:  -3 H:  -4 V: -11>
// <D: -13 H: -24 V: -14>                                                                      <D:  -8 H: -18 V: -16>   // <D:  -8 H: -18 V: -16>                                                                      <D:   1 H:   0 V: -10>   // <D:   1 H:   0 V: -10>                                                                      <D:   2 H:   0 V:  -6>  // <D:   2 H:   0 V:  -6>                                                                      <D:   1 H:   0 V: -10>
// <D: -14 H: -25 V: -15> <D:  -9 H: -20 V: -10> <D: -13 H: -21 V: -14> <D: -15 H: -22 V: -16> <D: -11 H: -22 V: -17>   // <D:   3 H:  -8 V:  -3> <D:   3 H:  -8 V:   2> <D:   9 H:  -2 V:  -2> <D:   8 H:  -3 V:  -3> <D:   6 H:  -4 V:  -5>   // <D:  10 H:   0 V:  -1> <D:   0 H:  -1 V:  -8> <D:  -1 H:  -2 V:  -9> <D:  -2 H:  -3 V:  -6> <D:  -3 H:  -4 V:  -7>  // <D:   1 H:   0 V:  -3> <D:   3 H:  -1 V:   2> <D:   9 H:  -2 V:  -2> <D:   8 H:  -3 V:  -3> <D:   6 H:  -4 V:  -5>

// +14                                                                                                                  // -3                                                                                                                   // -10                                                                                                                 // -1
// <D: -14 H: -25 V: -15> <D:  -9 H: -20 V: -10> <D: -13 H: -21 V: -14> <D: -15 H: -22 V: -16> <D: -11 H: -22 V: -17>   // <D:   3 H:  -8 V:  -3> <D:   3 H:  -8 V:   2> <D:   9 H:  -2 V:  -2> <D:   8 H:  -3 V:  -3> <D:   3 H:  -7 V:  -8>   // <D:  10 H:   0 V:  -1> <D:   0 H:  -1 V:  -8> <D:  -1 H:  -2 V:  -9> <D:  -2 H:  -3 V:  -6> <D: -13 H: -14 V: -17>  // <D:   0 H:  -1 V:  -4> <D:   2 H:  -2 V:   1> <D:   8 H:  -3 V:  -3> <D:   7 H:  -4 V:  -4> <D:   5 H:  -5 V:  -6>
// <D: -15 H: -26 V: -16>                                                                      <D: -15 H: -16 V: -18>   // <D:  -3 H:  -4 V:  -4>                                                                      <D:   0 H: -11 V:  -9>   // <D:   0 H: -11 V:  -9>                                                                      <D:  -9 H: -10 V: -18>  // <D:   4 H:   3 V:  -5> <D:   3 H:   2 V:   0> <D:   2 H:   1 V:  -4> <D:   3 H:   0 V:  -5> <D:   2 H:  -1 V:  -7>
// <D: -16 H: -27 V: -17>                                                                      <D: -10 H: -13 V: -19>   // <D:   2 H:   1 V:  -5>                                                                      <D:  -4 H:  -6 V: -10>   // <D:  -4 H:  -6 V: -10>                                                                      <D:  -4 H: -13 V: -15>  // <D:   9 H:   0 V:  -2> <D:   0 H:  -1 V:  -1> <D:  -1 H:  -2 V:  -5> <D:  -2 H:  -3 V:  -6> <D:  -2 H:  -4 V:  -8>
// <D: -17 H: -28 V: -18>                                                                      <D:  -6 H: -17 V: -17>   // <D:   8 H:  -3 V:  -3>                                                                      <D:  -9 H: -10 V: -11>   // <D:  -9 H: -10 V: -11>                                                                      <D:  -7 H: -17 V: -16>  // <D:   6 H:  -4 V:  -3> <D:  13 H:   2 V:   2> <D:   4 H:  -3 V:  -6> <D:  -3 H:  -4 V:  -7> <D:   0 H:  -1 V:  -9>
// <D: -18 H: -29 V: -19> <D:  -8 H: -16 V: -14> <D: -12 H: -17 V: -18> <D: -13 H: -18 V: -14> <D: -17 H: -19 V: -18>   // <D:  -6 H: -14 V:  -7> <D:   0 H: -11 V:  -5> <D:  -2 H: -12 V:  -9> <D:   6 H:  -5 V:  -5> <D:  -4 H:  -6 V: -12>   // <D:  -7 H:  -9 V: -15> <D:  -9 H: -10 V: -10> <D: -10 H: -11 V: -14> <D: -11 H: -12 V: -15> <D: -11 H: -13 V: -17>  // <D:   2 H:   0 V:  -4> <D:   2 H:  -1 V:   1> <D:   8 H:  -2 V:  -3> <D:   8 H:  -3 V:  -3> <D:   5 H:  -4 V:  -6>

// +4                                                                                                                   // +6                                                                                                                   // +7                                                                                                                  // -2
// <D:  -4 H: -15 V:  -5> <D:   1 H: -10 V:   0> <D:  -3 H: -11 V:  -4> <D:   1 H: -10 V:   0> <D:  -3 H: -11 V:  -4>   // <D:  -6 H: -14 V:  -7> <D:   0 H: -11 V:  -5> <D:  -2 H: -12 V:  -9> <D:   6 H:  -5 V:  -5> <D:   2 H:   0 V:  -6>   // <D:  -7 H:  -9 V: -15> <D:  -9 H: -10 V: -10> <D: -10 H: -11 V: -14> <D: -11 H: -12 V: -15> <D:  -4 H:  -6 V: -10>  // <D:   2 H:   0 V:  -4> <D:   2 H:  -1 V:   1> <D:   8 H:  -2 V:  -3> <D:   8 H:  -3 V:  -3> <D:   3 H:  -6 V:  -8>
// <D:  -1 H: -12 V:  -2> <D:   4 H:  -7 V:   3> <D:   9 H:  -2 V:  -1> <D:   4 H:  -3 V:   3> <D:   0 H:  -4 V:  -1>   // <D:  -1 H:  -5 V:  -2> <D:   1 H:  -6 V:   0> <D:   1 H:  -7 V:  -4> <D:   1 H:  -8 V:   0> <D:   7 H:  -4 V:  -4>   // <D:   5 H:  -6 V:  -6> <D:   4 H:  -7 V:  -4> <D:   2 H:  -8 V:  -8> <D:  -8 H:  -9 V:  -9> <D:  -9 H: -10 V: -11>  // <D:  -5 H:  -6 V:  -7> <D:  -1 H:  -7 V:  -2> <D:  -5 H:  -8 V:  -6> <D:   1 H:  -9 V:  -6> <D:   1 H: -10 V:  -9>
//
