// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_with_trace_recursion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the alignment algorithm using affine gap costs with trace
 *        information.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 *
 * \todo Maybe document trace_t. But we might also get rid of the template argument completely. Check this later.
 */
template <typename score_t, typename trace_t, bool truncate_score>
//!\cond
    requires (arithmetic<score_t> || simd_concept<score_t>) &&
             (std::same_as<trace_t, trace_directions> || simd_concept<trace_t>)
//!\endcond
class policy_affine_gap_with_trace_recursion : protected policy_affine_gap_recursion<score_t, truncate_score>
{
public:
    //!\brief The type of the base policy.
    using base_t = policy_affine_gap_recursion<score_t, truncate_score>;

    // Import base types.
    using typename base_t::score_type;
    using typename base_t::affine_score_tuple_t;

    //!\brief The trace type to use.
    using trace_type = trace_t;
    //!\brief The internal tuple storing the trace directions of an affine cell.
    using affine_trace_tuple_t = std::tuple<trace_type, trace_type, trace_type>;
    //!\brief The affine cell type returned by the functions.
    using affine_cell_type = affine_cell_proxy<std::pair<affine_score_tuple_t, affine_trace_tuple_t>>;

    // Import base member.
    using base_t::gap_extension_score;
    using base_t::gap_open_score;
    using base_t::first_row_is_free;
    using base_t::first_column_is_free;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_with_trace_recursion() = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion(policy_affine_gap_with_trace_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion(policy_affine_gap_with_trace_recursion &&) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion & operator=(policy_affine_gap_with_trace_recursion const &)
        = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion & operator=(policy_affine_gap_with_trace_recursion &&)
        = default; //!< Defaulted.
    ~policy_affine_gap_with_trace_recursion() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    explicit policy_affine_gap_with_trace_recursion(alignment_configuration_t const & config) : base_t{config}
    {}
    //!\}

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::compute_inner_cell
    template <typename affine_cell_t>
    affine_cell_type compute_inner_cell(score_type diagonal_score,
                                        affine_cell_t previous_cell,
                                        score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();
        trace_type best_trace = this->maybe_convert_to_simd(trace_directions::diagonal);

        auto cmp_mask = (diagonal_score < vertical_score);
        diagonal_score = cmp_mask ? vertical_score : diagonal_score;
        best_trace = cmp_mask ? previous_cell.vertical_trace() : (best_trace | previous_cell.vertical_trace());

        cmp_mask = (diagonal_score < horizontal_score);
        diagonal_score = cmp_mask ? horizontal_score : diagonal_score;
        best_trace = cmp_mask ? previous_cell.horizontal_trace() : (best_trace | previous_cell.horizontal_trace());

        this->truncate_score_below_zero(diagonal_score, best_trace);

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        cmp_mask = (vertical_score < tmp);
        vertical_score = cmp_mask ? tmp : vertical_score;
        trace_type next_vertical_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::up_open)
                                                  : this->maybe_convert_to_simd(trace_directions::up);

        cmp_mask = (horizontal_score < tmp);
        horizontal_score = cmp_mask ? tmp : horizontal_score;
        trace_type next_horizontal_trace = cmp_mask ? this->maybe_convert_to_simd(trace_directions::left_open)
                                                    : this->maybe_convert_to_simd(trace_directions::left);

        return {{diagonal_score, horizontal_score, vertical_score},
                {best_trace, next_horizontal_trace, next_vertical_trace}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_origin_cell
    affine_cell_type initialise_origin_cell() const noexcept
    {
        return {base_t::initialise_origin_cell(),
                {this->maybe_convert_to_simd(trace_directions::none),
                 first_row_is_free ? this->maybe_convert_to_simd(trace_directions::none)
                                   : this->maybe_convert_to_simd(trace_directions::left_open),
                 first_column_is_free ? this->maybe_convert_to_simd(trace_directions::none)
                                      : this->maybe_convert_to_simd(trace_directions::up_open)}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_column_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        return {base_t::initialise_first_column_cell(previous_cell),
                {previous_cell.vertical_trace(),
                 this->maybe_convert_to_simd(trace_directions::left_open),
                 first_column_is_free ? this->maybe_convert_to_simd(trace_directions::none)
                                      : this->maybe_convert_to_simd(trace_directions::up)}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_row_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        return {base_t::initialise_first_row_cell(previous_cell),
                {previous_cell.horizontal_trace(),
                 first_row_is_free ? this->maybe_convert_to_simd(trace_directions::none)
                                   : this->maybe_convert_to_simd(trace_directions::left),
                 this->maybe_convert_to_simd(trace_directions::up_open)}};
    }
};
} // namespace seqan3::detail
