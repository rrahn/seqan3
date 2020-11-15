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

#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_simd_adapter_iterator.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/simd/views/to_simd.hpp>

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
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
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

            if constexpr (traits_type::is_debug)
                this->initialise_debug_matrices(sequence1_size, sequence2_size);

            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair), alignment_matrix, index_matrix);
            this->make_result_and_invoke(*this,
                                         std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         this->optimal_score,
                                         this->optimal_coordinate,
                                         alignment_builder(alignment_matrix),
                                         callback);
        }
    }
    //!\}

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires traits_type::is_vectorised && std::invocable<callback_t, alignment_result_type>
    //!\endcond
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        auto && [simd_seq1_collection, simd_seq2_collection] = preprocess_sequences(indexed_sequence_pairs);

        auto && [alignment_matrix, index_matrix] = this->acquire_matrices(simd_seq1_collection.size(),
                                                                          simd_seq2_collection.size());

        compute_matrix(simd_seq1_collection, simd_seq2_collection, alignment_matrix, index_matrix);

        postprocess_result(std::forward<indexed_sequence_pairs_t>(indexed_sequence_pairs),
                           alignment_matrix,
                           std::forward<callback_t>(callback));
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
        static_assert(traits_type::is_vectorised, "This preprocess function is only callable in vectorisation mode.");

        using simd_collection_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;

        // Extract the batch of sequences for the first and the second sequence.
        auto seq1_collection = indexed_sequence_pairs | views::get<0> | views::get<0>;
        auto seq2_collection = indexed_sequence_pairs | views::get<0> | views::get<1>;

        this->initialise_tracker(seq1_collection, seq2_collection);

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
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename alignment_matrix_t, typename callback_t>
    void postprocess_result(indexed_sequence_pairs_t && indexed_sequence_pairs,
                            alignment_matrix_t && alignment_matrix,
                            callback_t && callback)
    {
        static_assert(traits_type::is_vectorised, "This preprocess function is only callable in vectorisation mode.");
        // Check that the sequence collection and the length of the vector have indeed the same size.
        assert(std::ranges::distance(indexed_sequence_pairs) <= simd::simd_traits<score_type>::length);

        using original_score_t = typename traits_type::original_score_type;

        size_t index = 0;
        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            original_score_t const padding_offset = this->padding_offsets[index];
            original_score_t score = this->optimal_score[index] -
                                     (padding_offset * this->scoring_scheme.padding_match_score());
            matrix_coordinate coordinate{
                row_index_type{static_cast<size_t>(this->optimal_coordinate.row[index] - padding_offset)},
                column_index_type{static_cast<size_t>(this->optimal_coordinate.col[index] - padding_offset)}
            };

            this->make_result_and_invoke(*this,
                                         std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         std::move(score),
                                         std::move(coordinate),
                                         alignment_builder(alignment_matrix, index),
                                         callback);
            ++index;
        }
    }

    /*!\brief Returns an invocable to build the aligned sequence if requested.
     *
     * \tparam alignment_matrix_t The type of the alignment matrix used to generate an aligned sequence builder for.
     *
     * \param[in] alignment_matrix The alignment matrix.
     * \param[in] simd_lane The lane of the simd vector to generate the aligned sequences for.
     *
     * \returns Either an invocable around a seqan3::detail::aligned_sequence_builder.
     *
     * \details
     *
     * This function returns an invocable wrapping the call to the seqan3::detail::aligned_sequence_builder.
     * The returned invocable is used by the seqan3::detail::policy_alignment_result_builder to compute the begin
     * positions or the aligned sequences if requested by the alignment configuration.
     * The signature of the returned invocable looks like the following:
     *
     * ```cpp
     * fn(auto && sequence1, auto && sequence2, auto coordinate)
     * ```
     * The first two arguments are the sequences to build the aligned sequences for. The third parameter is the
     * coordinate pointing to the begin of the trace to follow, which is associated with the computed alignment optimum.
     *
     * ### Complexity
     *
     * Constant.
     *
     * Calling the returned invocable is linear in the length of the generated trace: \f$ O(n+m) \f$, with
     * \f$ n,m \f$ being the size of the first, respectively second sequence.
     *
     * ### Exception
     *
     * The returned invocable might throw the following exceptions:
     * * std::out_of_range error, if the given coordinate exceeds the dimensions of the underlying trace matrix.
     * * seqan3::invalid_alignment_configuration, if the invocable is called when the alignment was configured without
     *   computing additional trace information.
     */
    template <typename alignment_matrix_t>
    auto alignment_builder([[maybe_unused]] alignment_matrix_t const & alignment_matrix, size_t const simd_lane = 0)
    {
        // Return a callable that can be invoked by the sequence builder.
        return [&, simd_lane] (auto && sequence1, auto && sequence2, auto coordinate)
        {
            // clip upper_diagonal
            int32_t upper_diagonal = std::clamp<int32_t>(this->upper_diagonal, 0, std::ranges::distance(sequence1));

            // If banded alignment, the row coordinate is mapped from the global matrix coordinate to the
            // internally defined one.
            if constexpr (traits_type::is_banded)
            {
                int32_t column_coordinate = static_cast<int32_t>(coordinate.col);
                coordinate.row -= (column_coordinate > upper_diagonal) * (column_coordinate - upper_diagonal);
            }

            if constexpr (traits_type::requires_trace_information)
            {
                // Get the matrix iterator at the specified alignment coordinate.
                auto trace_matrix_iter = alignment_matrix.matrix_iterator_at(coordinate);

                using matrix_iter_t = decltype(trace_matrix_iter);
                using matrix_iter_adapter_t =
                    lazy_conditional_t<traits_type::is_vectorised,
                                       lazy<trace_matrix_simd_adapter_iterator, matrix_iter_t>,
                                       matrix_iter_t>;
                using trace_iterator_t = trace_iterator<matrix_iter_adapter_t>;
                using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

                // Create the builder and return the generated alignment.
                aligned_sequence_builder builder{std::forward<decltype(sequence1)>(sequence1),
                                                 std::forward<decltype(sequence2)>(sequence2)};

                auto get_adapter_iter = [&] (auto iter)
                {
                    if constexpr (traits_type::is_vectorised)
                        return matrix_iter_adapter_t{std::move(iter), simd_lane};
                    else
                        return iter;
                };

                return builder(path_t{trace_iterator_t{get_adapter_iter(std::move(trace_matrix_iter)),
                                                       column_index_type{upper_diagonal}},
                                      std::default_sentinel});
            }
            else
            {
                throw seqan3::invalid_alignment_configuration{"You are trying to invoke the aligned sequence builder, "
                                                              "but the selected configuration disables the computation "
                                                              "of the trace."};
            }
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
        assert(static_cast<size_t>(std::ranges::distance(sequences)) <= traits_type::alignments_per_vector);

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

        if constexpr (traits_type::is_debug)
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
        using score_type = typename traits_type::score_type;

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

        if constexpr (traits_type::is_debug)
            this->log_alignment_matrix_column(cell_index_column, alignment_column
                                                               | views::take(std::ranges::distance(sequence2) + 1));
    }
};

} // namespace seqan3::detail
