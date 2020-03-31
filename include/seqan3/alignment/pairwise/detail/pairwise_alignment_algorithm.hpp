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

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

template <typename alignment_config_t, typename ...policies_t>
class pairwise_alignment_algorithm : protected policies_t...
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using scoring_scheme_type =  typename traits_type::scoring_scheme_t;

    scoring_scheme_type scoring_scheme;

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

    pairwise_alignment_algorithm(alignment_config_t const & config) : policies_t{config}...
    {
        scoring_scheme = seqan3::get<align_cfg::scoring>(config).value;
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given range over indexed sequence pairs.
     * \tparam indexed_sequence_pairs_t The type of indexed_sequence_pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     *
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     *
     * \returns A std::vector over seqan3::alignment_result with the requested alignment results for every
     *          sequence pair in the given range.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment for each
     * sequence pair. The space and runtime complexities depend on the selected configurations (see below).
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
    template <indexed_sequence_pair_range indexed_sequence_pairs_t>
    //!\cond
        requires !traits_type::is_vectorised
    //!\endcond
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs)
    {
        using sequence_pairs_t = std::tuple_element_t<0, std::ranges::range_value_t<indexed_sequence_pairs_t>>;
        using sequence1_t = std::tuple_element_t<0, sequence_pairs_t>;
        using sequence2_t = std::tuple_element_t<1, sequence_pairs_t>;
        using result_t = typename align_result_selector<sequence1_t, sequence2_t, alignment_config_t>::type;

        using std::get;

        std::vector<alignment_result<result_t>> results{};
        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            result_t res{};
            res.id = idx;
            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));

            if constexpr (traits_type::result_type_rank >= 0)
                res.score = this->optimum_score();

            if constexpr (traits_type::result_type_rank >= 1)
                res.back_coordinate = alignment_coordinate{column_index_type{this->optimum_coordinate().column_index()},
                                                           row_index_type{this->optimum_coordinate().row_index()}};

            results.emplace_back(std::move(res));
        }
        return results;
    }

protected:

    /*!\brief Compute the alignment by iterating over the alignment matrix in a column wise manner.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        this->reset_tracker();
        auto alignment_matrix = this->acquire_alignment_matrix(sequence1, sequence2);
        auto matrix_iterator = alignment_matrix.begin();

        initialise_column(matrix_iterator);

        // ---------------------------------------------------------------------
        // Recursion phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        auto track_cell = [this] (auto && ...args)
        {
            return this->track_cell(std::forward<decltype(args)>(args)...);
        };

        auto track_last_row_cell = [this] (auto && ...args)
        {
            return this->track_last_row_cell(std::forward<decltype(args)>(args)...);
        };

        ++matrix_iterator; // Move to next column.
        auto seq1_iterator = sequence1.begin();
        const size_t seq1_size = std::ranges::distance(sequence1);

        for (size_t column_index = 1; column_index < seq1_size; ++column_index, ++seq1_iterator, ++matrix_iterator)
            compute_column(matrix_iterator, seq1_iterator, sequence2.begin(), track_cell, track_last_row_cell);

        // ---------------------------------------------------------------------
        // Final phase: compute last column and report optimum.
        // ---------------------------------------------------------------------

        auto track_last_column_cell = [this] (auto && ...args)
        {
            return this->track_last_column_cell(std::forward<decltype(args)>(args)...);
        };

        auto track_final_cell = [this] (auto && ...args)
        {
            return this->track_final_cell(std::forward<decltype(args)>(args)...);
        };

        compute_column(matrix_iterator, seq1_iterator, sequence2.begin(), track_last_column_cell, track_final_cell);
    }

    template <std::forward_iterator matrix_iter_t>
    //!\cond
        requires std::ranges::sized_range<std::tuple_element_t<0, std::iter_reference_t<matrix_iter_t>>>
    //!\endcond
    void initialise_column(matrix_iter_t matrix_iter)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto [first_column, coordinate_column] = *matrix_iter;
        auto first_column_iter = std::ranges::begin(first_column);
        auto coordinate_iter = std::ranges::begin(coordinate_column);

        *first_column_iter = this->track_cell(this->initialise_origin_cell(), *coordinate_iter);

        // ---------------------------------------------------------------------
        // Recursion phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        ++first_column_iter;
        ++coordinate_iter;
        for (size_t i = 1; i < std::ranges::size(first_column) - 1; ++i, ++first_column_iter, ++coordinate_iter)
            *first_column_iter = this->track_cell(this->initialise_first_column_cell(*first_column_iter),
                                                  *coordinate_iter);

        // ---------------------------------------------------------------------
        // Final phase: compute last cell in column
        // ---------------------------------------------------------------------

        *first_column_iter = this->track_last_row_cell(this->initialise_first_column_cell(*first_column_iter),
                                                       *coordinate_iter);
    }

    template <std::forward_iterator matrix_iter_t,
              std::readable sequence1_iterator_t,
              std::forward_iterator sequence2_iterator_t,
              std::invocable<
                std::ranges::range_value_t<
                    std::tuple_element_t<0, std::iter_reference_t<matrix_iter_t>>>,
                std::ranges::range_reference_t<
                    std::tuple_element_t<1, std::iter_reference_t<matrix_iter_t>>>> track_cell_fn_t,
              std::invocable<
                std::ranges::range_value_t<
                    std::tuple_element_t<0, std::iter_reference_t<matrix_iter_t>>>,
                std::ranges::range_reference_t<
                    std::tuple_element_t<1, std::iter_reference_t<matrix_iter_t>>>> track_last_row_cell_fn_t>
    //!\cond
        requires std::ranges::sized_range<std::tuple_element_t<0, std::iter_reference_t<matrix_iter_t>>>
    //!\endcond
    void compute_column(matrix_iter_t matrix_iter,
                        sequence1_iterator_t seq1_iter,
                        sequence2_iterator_t seq2_iter,
                        track_cell_fn_t && track_cell,
                        track_last_row_cell_fn_t && track_last_row_cell)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------
        auto [alignment_column, coordinate_column] = *matrix_iter;
        auto column_iter = std::ranges::begin(alignment_column);
        auto coordinate_iter = std::ranges::begin(coordinate_column);

        auto cell = *column_iter;
        typename traits_type::score_t diagonal = cell.optimal_score();
        *column_iter = track_cell(this->initialise_first_row_cell(cell), *coordinate_iter);

        // ---------------------------------------------------------------------
        // Recursion phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        ++column_iter;
        ++coordinate_iter;
        const size_t column_size = std::ranges::size(alignment_column) - 1;
        for (size_t row_index = 1; row_index < column_size; ++row_index, ++seq2_iter, ++column_iter, ++coordinate_iter)
        {
            auto cell = *column_iter;
            typename traits_type::score_t next_diagonal = cell.optimal_score();
            *column_iter = track_cell(this->compute_inner_cell(diagonal,
                                                               cell,
                                                               this->scoring_scheme.score(*seq1_iter, *seq2_iter)),
                                      *coordinate_iter);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: compute last cell in column
        // ---------------------------------------------------------------------

        *column_iter = track_last_row_cell(this->compute_inner_cell(diagonal,
                                                                    *column_iter,
                                                                    this->scoring_scheme.score(*seq1_iter, *seq2_iter)),
                                           *coordinate_iter);
    }
};
} // namespace seqan3::detail
