// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_kernel.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/simd/concept.hpp>

namespace seqan3::detail
{

template <typename recursion_policy_t, typename scoring_policy_t, typename tracker_policy_t>
struct pairwise_alignment_state : public recursion_policy_t, public scoring_policy_t, public tracker_policy_t
{
public:
    pairwise_alignment_state() = default;
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    pairwise_alignment_state(alignment_configuration_t const & config) :
        recursion_policy_t{config},
        scoring_policy_t{config},
        tracker_policy_t{config}
    {}
};

// When are they initialised?
// template <typename recursion_policy_t, typename scoring_policy_t, typename tracker_policy_t>
struct pairwise_alignment_kernel //: public recursion_policy_t, public scoring_policy_t, public tracker_policy_t
{
public:

    pairwise_alignment_kernel() = default;
    // template <typename alignment_configuration_t>
    // //!\cond
    //     requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    // //!\endcond
    // pairwise_alignment_kernel(alignment_configuration_t const & config) :
    //     recursion_policy_t{config},
    //     scoring_policy_t{config},
    //     tracker_policy_t{config}
    // {}

    // how to communicate with a logger from outside?
    // Don't like to provide more data
    // We should be able to ask

    template <typename sequence1_t, typename sequence2_t, typename alignment_matrix_t, typename index_matrix_t, typename state_t>
    constexpr auto operator()(sequence1_t && sequence1,
                              sequence2_t && sequence2,
                              alignment_matrix_t & alignment_matrix,
                              index_matrix_t && index_matrix,
                              state_t & state) const // are we sure we do not modify any state?
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        state.reset_optimum(); // Reset the tracker for the new alignment computation.

        auto alignment_matrix_it = alignment_matrix.begin();
        auto indexed_matrix_it = index_matrix.begin();

        initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2, state);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        for (auto alphabet1 : sequence1)
            compute_column(*++alignment_matrix_it,
                           *++indexed_matrix_it,
                           state.scoring_scheme_profile_column(alphabet1),
                           sequence2,
                           state);

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        auto && alignment_column = *alignment_matrix_it;
        auto && cell_index_column = *indexed_matrix_it;

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        // TODO: not working in const version. -> so we either need non-const
        // or can we pass this as a state into the algorithm?
        state.track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        for ([[maybe_unused]] auto && unused : sequence2)
            state.track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        state.track_final_cell(*alignment_column_it, *cell_index_column_it);

        // return state;
    }

protected:
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
              std::ranges::forward_range sequence2_t,
              typename state_t>
    void initialise_column(alignment_column_t && alignment_column,
                           cell_index_column_t && cell_index_column,
                           sequence2_t && sequence2,
                           state_t & state) const
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto first_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();
        *first_column_it = state.track_cell(state.initialise_origin_cell(), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for ([[maybe_unused]] auto const & unused : sequence2)
        {
            ++first_column_it;
            *first_column_it = state.track_cell(state.initialise_first_column_cell(*first_column_it),
                                                *cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell of initial column
        // ---------------------------------------------------------------------

        state.track_last_row_cell(*first_column_it, *cell_index_column_it);

        // TODO: enable
        // this->log_alignment_matrix_column(cell_index_column, alignment_column
        //                                                    | views::take(std::ranges::distance(sequence2)+ 1));
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
              std::ranges::forward_range sequence2_t,
              typename state_t>
    //!\cond
        requires semialphabet<alphabet1_t> || simd_concept<alphabet1_t>
    //!\endcond
    void compute_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        alphabet1_t const & alphabet1,
                        sequence2_t && sequence2,
                        state_t & state) const
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        auto diagonal = cell.best_score();
        *alignment_column_it = state.track_cell(state.initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto const & alphabet2 : sequence2)
        {
            auto cell = *++alignment_column_it;
            auto next_diagonal = cell.best_score();
            *alignment_column_it = state.track_cell(
                state.compute_inner_cell(diagonal, cell, state.scoring_scheme.score(alphabet1, alphabet2)),
                *cell_index_column_it);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        state.track_last_row_cell(*alignment_column_it, *cell_index_column_it);

        // this->log_alignment_matrix_column(cell_index_column, alignment_column
        //                                                    | views::take(std::ranges::distance(sequence2) + 1));
    }
};


}  // namespace seqan3::detail
