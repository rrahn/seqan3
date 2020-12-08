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

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 *
 * \tparam score_t The configured score type.
 * \tparam result_t The configured result type.
 * \tparam policies_t Variadic template argument for the different policies of this alignment algorithm.
 *
 * \details
 *
 * The algorithm computes a classical dynamic programming matrix to obtain the alignment for a pair of sequences.
 * The matrix is computed column-wise and expects that all policies adhere to this computation layout.
 * The algorithm can also compute alignments for batches of seuence pairs concurrently using vectorisation on
 * dedicated SIMD registers.
 * After the computation a user defined callback function is invoked for every generated seqan3::alignment_result.
 * This means that one seqan3::alignment_result is associated with only one result, even if a batch of alignments is
 * computed in the vetorised mode.
 */
template <typename score_t, typename result_t,
          typename state_t,
          typename ...policies_t>
class pairwise_alignment_algorithm : protected policies_t...
{
protected:

    // ----------------------------------------------------------------------------
    // static member
    // ----------------------------------------------------------------------------

    //!\brief Flag indicating whether the alignment is executed in vectorised mode.
    static constexpr bool is_vectorised = simd_concept<score_t>;

    // ----------------------------------------------------------------------------
    // type definition
    // ----------------------------------------------------------------------------

    //!\brief The configured score type.
    using score_type = score_t;
    //!\brief Helper template alias to extract the scalar type of a simd type in a lazy conditional.
    template <simd_concept simd_score_t>
    using scalar_type_of = typename simd_traits<simd_score_t>::scalar_type;
    //!\brief The original score type.
    using original_score_type = lazy_conditional_t<is_vectorised, lazy<scalar_type_of, score_type>, score_type>;
    //!\brief The configured alignment result type.
    using alignment_result_type = result_t;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");

    state_t alignment_state{};

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
     *
     * \tparam alignment_configuration_t The type of the alignment configuration.
     *
     * \param[in] config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the base policies of the alignment algorithm.
     */
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    pairwise_alignment_algorithm(alignment_configuration_t const & config) :
        policies_t(config)...,
        alignment_state{config}
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
     * on the given \ref seqan3_align_cfg_output_configurations "seqan3::align_cfg::output_*" per sequence pair.
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
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            size_t const sequence1_size = std::ranges::distance(get<0>(sequence_pair));
            size_t const sequence2_size = std::ranges::distance(get<1>(sequence_pair));

            auto && [alignment_matrix, index_matrix] = this->acquire_matrices(sequence1_size, sequence2_size);

            this->initialise_debug_matrices(sequence1_size, sequence2_size);

            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair), alignment_matrix, index_matrix);
            this->make_result_and_invoke(*this,
                                         std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         this->optimal_score,
                                         this->optimal_coordinate,
                                         trace_path_generator(alignment_matrix),
                                         callback);
        }
    }
    //!\}

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires is_vectorised && std::invocable<callback_t, alignment_result_type>
    //!\endcond
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        auto && [simd_seq1_collection, simd_seq2_collection] = preprocess_sequences(indexed_sequence_pairs);

        auto && [alignment_matrix, index_matrix] = this->acquire_matrices(simd_seq1_collection.size(),
                                                                          simd_seq2_collection.size());
        // compute_matrix(simd_seq1_collection, simd_seq2_collection, alignment_matrix, index_matrix);
        // pairwise_alignment_kernel{}{alignment_kernel};
        pairwise_alignment_kernel{}(simd_seq1_collection,
                                    simd_seq2_collection,
                                    alignment_matrix,
                                    index_matrix,
                                    alignment_state);

        postprocess_result(std::forward<indexed_sequence_pairs_t>(indexed_sequence_pairs),
                           alignment_matrix,
                           std::forward<callback_t>(callback),
                           alignment_state);
    }

protected:
    /*!\brief Preprocess to transform the sequence collections to simd vectors.
     * \tparam indexed_sequence_pairs_t The type of the simd sequence pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     *
     * \returns A tuple over the first and second sequence collection transformed to vectors of simd types.
     *
     * \details
     *
     * Decomposes the indexed sequence pairs into a collection for the first sequences and the second sequences.
     * Each collection is converted to a vector over simd data types according to configurations set in the
     * initialisation of the algorithm. Both vectors are thread local variables and are cached between
     * invocations of the alignment algorithm in order to reduce memory allocations.
     * Finally, a tuple with references to the thread local simd vectors is returned.
     *
     * \throws std::bad_alloc if the allocation of the simd vectors exceeds the available memory.
     */
    template <indexed_sequence_pair_range indexed_sequence_pairs_t>
    auto preprocess_sequences(indexed_sequence_pairs_t && indexed_sequence_pairs)
    {
        static_assert(is_vectorised, "This preprocess function is only callable in vectorisation mode.");

        using simd_collection_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;

        // Extract the batch of sequences for the first and the second sequence.
        auto seq1_collection = indexed_sequence_pairs | views::get<0> | views::get<0>;
        auto seq2_collection = indexed_sequence_pairs | views::get<0> | views::get<1>;

        // this->initialise_tracker(seq1_collection, seq2_collection);
        alignment_state.initialise_tracker(seq1_collection, seq2_collection);

        // Convert batch of sequences to sequence of simd vectors.
        thread_local simd_collection_t simd_seq1_collection{};
        thread_local simd_collection_t simd_seq2_collection{};

        convert_batch_of_sequences_to_simd_vector(simd_seq1_collection,
                                                  seq1_collection,
                                                  this->scoring_scheme.padding_symbol);
        convert_batch_of_sequences_to_simd_vector(simd_seq2_collection,
                                                  seq2_collection,
                                                  this->scoring_scheme.padding_symbol);
        return std::tie(simd_seq1_collection, simd_seq2_collection);
    }

    /*!\brief Postprocesses the result computation of the vectorised alignment.
     * \tparam indexed_sequence_pairs_t The type of the simd sequence pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \tparam alignment_matrix_t The type of the alignment matrix created for the this alignment algorithm.
     * \tparam callback_t The type of the callback to be invoked on every alignment result.
     *
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     * \param[in] alignment_matrix The alignment matrix which is passed to the result builder.
     * \param[in] callback The callback which is invoked on the constructed alignment result.
     *
     * \details
     *
     * Iterates over the sequence pairs and extracts the corresponding score and coordinate for each alignment pair in
     * the collection. Note the process assumes that the size of the collection is less or equal than the size of the
     * simd vector. No further runtime checks are enforced and if the pre-condition is not met the program will be left
     * in the state of undefined behaviour.
     * Subsequently, the alignment result builder is invoked with the respective values to construct the alignment
     * result and invoke the callable with the new result.
     */
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename alignment_matrix_t, typename callback_t, typename tracker_t>
    void postprocess_result(indexed_sequence_pairs_t && indexed_sequence_pairs,
                            alignment_matrix_t && alignment_matrix,
                            callback_t && callback,
                            tracker_t const & tracker)
    {
        static_assert(is_vectorised, "This preprocess function is only callable in vectorisation mode.");
        // Check that the sequence collection and the length of the vector have indeed the same size.
        assert(std::ranges::distance(indexed_sequence_pairs) <= simd::simd_traits<score_type>::length);

        size_t index = 0;
        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            original_score_type const padding_offset = tracker.padding_offsets[index];
            original_score_type score = tracker.optimal_score[index] -
                                     (padding_offset * this->scoring_scheme.padding_match_score());
            matrix_coordinate coordinate{
                row_index_type{static_cast<size_t>(tracker.optimal_coordinate.row[index] - padding_offset)},
                column_index_type{static_cast<size_t>(tracker.optimal_coordinate.col[index] - padding_offset)}
            };

            this->make_result_and_invoke(*this,
                                         std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         std::move(score),
                                         std::move(coordinate),
                                         trace_path_generator(alignment_matrix, index),
                                         callback);
            ++index;
        }
    }

    /*!\brief Returns an invocable which generates a trace path over the given alignment matrix.
     *
     * \tparam alignment_matrix_t The type of the alignment matrix used to generate an aligned sequence builder for.
     *
     * \param[in] alignment_matrix The alignment matrix.
     * \param[in] simd_position The position of the simd vector to get the trace path for.
     *
     * \returns An invocable that, when invoked, returns the trace path starting at the given matrix coordinate.
     *
     * \details
     *
     * This function returns an invocable which returns the trace path from the trace matrix starting at the given
     * matrix coordinate. This trace path will be used later by the
     * seqan3::detail::policy_alignment_result_builder to generate the aligned sequences.
     * If the alignment was executed in vectorised mode, then the builder is invoked for each alignment stored in a
     * separate simd vector position. In the scalar mode this position is ignored.
     *
     * Following syntax can be expected from the returned invocable:
     * ```cpp
     * fn(auto coordinate) -> trace_path
     * ```
     * where the returned `trace_path` is a std::ranges::forward_range over seqan3::detail::trace_directions.
     *
     * ### Exception
     *
     * Might throw implementation defined errors.
     * If an error is thrown the state of this class is not changed (strong exception guarantee).
     */
    template <typename alignment_matrix_t>
    auto trace_path_generator(alignment_matrix_t const & alignment_matrix, size_t const simd_position = 0)
    {
        // Return a callable that can be invoked by the sequence builder.
        return [&, simd_position] (auto const & coordinate)
        {
            return this->trace_path_starting_at(alignment_matrix, coordinate, simd_position);
        };
    }

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
        assert(static_cast<size_t>(std::ranges::distance(sequences)) <= simd_traits<score_type>::length);

        simd_sequence.clear();
        for (auto && simd_vector_chunk : sequences | views::to_simd<score_type>(padding_symbol))
            std::ranges::move(simd_vector_chunk, std::cpp20::back_inserter(simd_sequence));
    }

    /*!\brief Compute the actual alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     * \tparam alignment_matrix_t The type of the alignment matrix; must model std::ranges::input_range and its
     *                            std::ranges::range_reference_t type must model std::ranges::forward_range.
     * \tparam index_matrix_t The type of the index matrix; must model std::ranges::input_range and its
     *                            std::ranges::range_reference_t type must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     * \param[in] alignment_matrix The alignment matrix to compute.
     * \param[in] index_matrix The index matrix corresponding to the alignment matrix.
     */
    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              std::ranges::input_range alignment_matrix_t,
              std::ranges::input_range index_matrix_t>
    //!\cond
        requires std::ranges::forward_range<std::ranges::range_reference_t<alignment_matrix_t>> &&
                 std::ranges::forward_range<std::ranges::range_reference_t<index_matrix_t>>
    //!\endcond
    void compute_matrix(sequence1_t && sequence1,
                        sequence2_t && sequence2,
                        alignment_matrix_t && alignment_matrix,
                        index_matrix_t && index_matrix)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        this->reset_optimum(); // Reset the tracker for the new alignment computation.

        auto alignment_matrix_it = alignment_matrix.begin();
        auto indexed_matrix_it = index_matrix.begin();

        initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        for (auto alphabet1 : sequence1)
            compute_column(*++alignment_matrix_it,
                           *++indexed_matrix_it,
                           this->scoring_scheme_profile_column(alphabet1),
                           sequence2);

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        auto && alignment_column = *alignment_matrix_it;
        auto && cell_index_column = *indexed_matrix_it;

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        this->track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        for ([[maybe_unused]] auto && unused : sequence2)
            this->track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        this->track_final_cell(*alignment_column_it, *cell_index_column_it);
    }

    /*!\brief Initialise the first column of the alignment matrix.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::forward_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
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
    template <std::ranges::forward_range alignment_column_t,
              std::ranges::forward_range cell_index_column_t,
              std::ranges::forward_range sequence2_t>
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

        this->log_alignment_matrix_column(cell_index_column, alignment_column
                                                           | views::take(std::ranges::distance(sequence2)+ 1));
    }

    /*!\brief Initialise any column of the alignment matrix except the first one.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::forward_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::forward_range.
     * \tparam alphabet1_t The type of the current symbol of sequence1.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] alphabet1 The current symbol of sequence1.
     * \param[in] sequence2 The second sequence to align against `alphabet1`.
     *
     * \details
     *
     * Computes the alignment for the given alignment matrix column. The function splits the computation of the column
     * into three phases: the initialisation phase, the iteration phase, and the final phase. In the initialisation
     * phase the first cell of the column is computed and in the iteration phase all remaining cells are computed.
     * In the final phase the last cell is possibly evaluated for a new alignment optimum.
     */
    template <std::ranges::forward_range alignment_column_t,
              std::ranges::forward_range cell_index_column_t,
              typename alphabet1_t,
              std::ranges::forward_range sequence2_t>
    //!\cond
        requires semialphabet<alphabet1_t> || simd_concept<alphabet1_t>
    //!\endcond
    void compute_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        alphabet1_t const & alphabet1,
                        sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        score_type diagonal = cell.best_score();
        *alignment_column_it = this->track_cell(this->initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto const & alphabet2 : sequence2)
        {
            auto cell = *++alignment_column_it;
            score_type next_diagonal = cell.best_score();
            *alignment_column_it = this->track_cell(
                this->compute_inner_cell(diagonal, cell, this->scoring_scheme.score(alphabet1, alphabet2)),
                *++cell_index_column_it);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*alignment_column_it, *cell_index_column_it);

        this->log_alignment_matrix_column(cell_index_column, alignment_column
                                                           | views::take(std::ranges::distance(sequence2) + 1));
    }
};

} // namespace seqan3::detail
