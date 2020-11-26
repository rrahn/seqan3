// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_with_trace_recursion_opt.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the alignment algorithm using affine gap costs with trace
 *        information.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename traits_type>
class policy_affine_gap_with_trace_recursion_opt : protected policy_affine_gap_recursion<traits_type>
{
protected:
    //!\brief The type of the base policy.
    using base_t = policy_affine_gap_recursion<traits_type>;

    // Import base types.
    using typename base_t::score_type;

    //!\brief The trace type to use.
    using trace_type = typename traits_type::trace_type;
    using affine_score_pair_type = std::pair<score_type, score_type>;
    using affine_score_cell_type = std::pair<affine_score_pair_type, affine_score_pair_type>;
    //!\brief The affine cell type returned by the functions.
    using affine_cell_type = affine_cell_proxy<std::pair<affine_score_cell_type, trace_type>>;

    // Import base member.
    using base_t::gap_extension_score;
    using base_t::gap_open_score;
    using base_t::first_row_is_free;
    using base_t::first_column_is_free;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_with_trace_recursion_opt() = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion_opt(policy_affine_gap_with_trace_recursion_opt const &) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion_opt(policy_affine_gap_with_trace_recursion_opt &&) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion_opt & operator=(policy_affine_gap_with_trace_recursion_opt const &)
        = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion_opt & operator=(policy_affine_gap_with_trace_recursion_opt &&)
        = default; //!< Defaulted.
    ~policy_affine_gap_with_trace_recursion_opt() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    explicit policy_affine_gap_with_trace_recursion_opt(alignment_configuration_t const & config) : base_t{config}
    {}
    //!\}

    template <typename affine_cell_t>
    affine_cell_type compute_inner_cell(score_type best_score, affine_cell_t previous_cell, score_type substitution_score)
    {
        //                best   , horizon ,  vertical            best_trace
        // The layout: {{{score_t, score_t}, {score_t, score_t}}, trace}
        auto [up_open, up_extension] = previous_cell.vertical_score();

        // Diagonal value + substitution score.
        best_score += substitution_score;

        // Compute the optimum coming from the upper cell
        up_open += gap_open_score;
        up_extension += gap_extension_score;

        auto cmp_mask = (up_open < up_extension);
        up_extension = cmp_mask ? up_extension : up_open;
        trace_type trace_up = cmp_mask ? this->maybe_convert_to_simd(trace_directions::up)
                                       : this->maybe_convert_to_simd(trace_directions::up_open);

        // Maximum between U and D.
        cmp_mask = (best_score < up_extension);
        best_score = cmp_mask ? up_extension : best_score;
        trace_type trace_best = cmp_mask ? trace_up
                                         : this->maybe_convert_to_simd(trace_directions::diagonal) | trace_up;

        // Compute the optimum coming from the left cell
        substitution_score = previous_cell.best_score() + gap_open_score;
        up_open = previous_cell.horizontal_score() + gap_extension_score;

        cmp_mask = (up_open < substitution_score);
        substitution_score = cmp_mask ? substitution_score : up_open;
        trace_up |= cmp_mask ? this->maybe_convert_to_simd(trace_directions::left)
                             : this->maybe_convert_to_simd(trace_directions::left_open);

        // Maximum between L and D.
        // Order of this has influence of precedence in trace back.
        // So here it is better to come from left than from up or diagonal.
        // This means if we get into this cell, we need to follow left instead
        // of up.
        // In a clean interface this would be coupled.
        // So the same structure that decides how to handle this must also
        // decide how to follow the trace. Or we need to document it carefully.
        // So the aligned sequence builder and this operation.
        cmp_mask = (best_score < substitution_score);
        best_score = cmp_mask ? substitution_score : best_score;
        trace_best = cmp_mask ? trace_up : trace_best | trace_up;

        // Now we have everything together!
        return {{{best_score, substitution_score}, {best_score, up_extension}}, trace_best};

        // 17 assignments
        // 5 additions
        // 4 comparisons
        // 8 blends
        // 2 or-operations
        // 36 operations
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::compute_inner_cell
    // template <typename affine_cell_t>
    // affine_cell_type compute_inner_cell(score_type diagonal_score,
    //                                     affine_cell_t previous_cell,
    //                                     score_type const sequence_score) const noexcept
    // {
    //     diagonal_score += sequence_score;
    //     score_type horizontal_score = previous_cell.horizontal_score();
    //     score_type vertical_score = previous_cell.vertical_score();
    //     trace_type best_trace = this->maybe_convert_to_simd(trace_directions::diagonal);

    //     auto cmp_mask = (diagonal_score < vertical_score);
    //     diagonal_score = cmp_mask ? vertical_score : diagonal_score;
    //     best_trace = cmp_mask ? previous_cell.vertical_trace() : (best_trace | previous_cell.vertical_trace());

    //     cmp_mask = (diagonal_score < horizontal_score);
    //     diagonal_score = cmp_mask ? horizontal_score : diagonal_score;
    //     best_trace = cmp_mask ? previous_cell.horizontal_trace() : (best_trace | previous_cell.horizontal_trace());

    //     this->truncate_score_below_zero(diagonal_score, best_trace);

    //     score_type tmp = diagonal_score + gap_open_score;
    //     vertical_score += gap_extension_score;
    //     horizontal_score += gap_extension_score;

    //     cmp_mask = (vertical_score < tmp);
    //     vertical_score = cmp_mask ? tmp : vertical_score;
    //     trace_type next_vertical_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::up_open)
    //                                               : this->maybe_convert_to_simd(trace_directions::up);

    //     cmp_mask = (horizontal_score < tmp);
    //     horizontal_score = cmp_mask ? tmp : horizontal_score;
    //     trace_type next_horizontal_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::left_open)
    //                                                 : this->maybe_convert_to_simd(trace_directions::left);

    //     return {{diagonal_score, horizontal_score, vertical_score},
    //             {best_trace, next_horizontal_trace, next_vertical_trace}};
    // }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_origin_cell
    affine_cell_type initialise_origin_cell() const noexcept
    {
        // {0, -inf, -inf}
        score_type horizontal_score = first_column_is_free ? score_type{} : gap_open_score;
        score_type vertical_score = first_row_is_free ? score_type{} : gap_open_score;
        return {{{score_type{}, horizontal_score}, {score_type{}, vertical_score}},
                this->maybe_convert_to_simd(trace_directions::none)};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_column_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        auto [up_open_score, up_extension_score] = previous_cell.vertical_score();
        up_open_score += gap_open_score;
        up_extension_score += gap_extension_score; // 0 + gap_open is better than gap_open + gap_extension as long as gap_extension < 0.

        auto cmp_mask = (up_extension_score < up_open_score);
        score_type best_score = cmp_mask ? up_open_score : up_extension_score;
        trace_type best_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::up_open)
                                         : this->maybe_convert_to_simd(trace_directions::up);
        best_score = first_column_is_free ? score_type{} : best_score;
        best_trace = first_column_is_free ? this->maybe_convert_to_simd(trace_directions::none) : best_trace;
        return {{{best_score, best_score + gap_open_score}, {best_score, best_score}}, best_trace};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_row_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        auto [left_open_score, left_extension_score] = previous_cell.vertical_score();
        left_open_score = previous_cell.best_score() + gap_open_score;
        left_extension_score = previous_cell.horizontal_score() + gap_extension_score; // 0 + gap_open is better than gap_open + gap_extension as long as gap_extension < 0.

        auto cmp_mask = (left_extension_score < left_open_score);
        score_type best_score = cmp_mask ? left_open_score : left_extension_score;
        trace_type best_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::left_open)
                                         : this->maybe_convert_to_simd(trace_directions::left);
        best_score = first_row_is_free ? score_type{} : best_score;
        best_trace = first_row_is_free ? this->maybe_convert_to_simd(trace_directions::none) : best_trace;
        return {{{best_score, best_score}, {best_score, best_score + gap_open_score}}, best_trace};

        // return {base_t::initialise_first_row_cell(previous_cell),
        //         {previous_cell.horizontal_trace(),
        //          first_row_is_free ? this->maybe_convert_to_simd(trace_directions::none)
        //                            : this->maybe_convert_to_simd(trace_directions::left),
        //          this->maybe_convert_to_simd(trace_directions::up_open)}};
    }
};
} // namespace seqan3::detail
