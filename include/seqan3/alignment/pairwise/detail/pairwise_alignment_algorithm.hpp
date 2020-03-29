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

#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/range/views/drop.hpp>

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
            res.score = compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));
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
    template <typename sequence1_t, typename sequence2_t>
    int32_t compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2)
    //!\cond
        requires !traits_type::is_banded
    //!\endcond
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        auto local_matrix = this->acquire_alignment_matrix(sequence1, sequence2);
        auto matrix_iterator = local_matrix.begin();
        auto alignment_column = *matrix_iterator;

        // Initialise the first column.
        *alignment_column.begin() = this->initialise_origin_cell();
        for (auto cell : alignment_column | seqan3::views::drop(1))
            cell = this->initialise_first_column_cell(cell);

        // ---------------------------------------------------------------------
        // Recursion phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        // Compute remaining matrix
        ++matrix_iterator; // Move to next column.
        for (auto it_col = sequence1.begin(); it_col != sequence1.end(); ++it_col, ++matrix_iterator)
        {
            alignment_column = *matrix_iterator;
            auto alignment_column_iter = alignment_column.begin();

            // Initialise first cell of current column.
            auto cell = *alignment_column_iter;
            typename traits_type::score_t diagonal = cell.optimal_score();
            *alignment_column_iter = this->initialise_first_row_cell(cell);

            ++alignment_column_iter;  // Move to next cell in column.
            for (auto it_row = sequence2.begin(); it_row != sequence2.end(); ++it_row, ++alignment_column_iter)
            {
                auto cell = *alignment_column_iter;
                typename traits_type::score_t next_diagonal = cell.optimal_score();
                *alignment_column_iter = this->compute_inner_cell(diagonal,
                                                                  cell,
                                                                  this->scoring_scheme.score(*it_col, *it_row));
                diagonal = next_diagonal;
            }
        }

        // ---------------------------------------------------------------------
        // Wrap up phase: track score of last column
        // ---------------------------------------------------------------------

        return alignment_column[std::ranges::distance(sequence2)].optimal_score();
    }
};
} // namespace seqan3::detail
