// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_recursion_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker_simd.hpp>

namespace seqan3::detail
{
/*!\brief Selects the recursion policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 *
 * \details
 *
 * Configures the recursion policy type based on the configured alignment algorithm traits.
 * The selector differentiates between banded alignments and alignments that need the trace information.
 */
template <typename algorithm_traits_t>
struct select_recursion_policy
{
private:
    //!\brief A flag indicating whether the trace must be computed.
    static constexpr bool with_trace = algorithm_traits_t::requires_trace_information;
    //!\brief A flag indicating whether the band must be computed.
    static constexpr bool with_band = algorithm_traits_t::is_banded;
    //!\brief A flag indicating whether the score must be truncated.
    static constexpr bool truncate_score = algorithm_traits_t::is_local;

    //!\brief The configured score type.
    using score_t = typename algorithm_traits_t::score_type;
    //!\brief The configured trace type.
    using trace_t = typename algorithm_traits_t::trace_type;

    //!\brief Selects the affine gap policy that is able to compute the trace.
    using policy_with_trace_t =
        std::conditional_t<with_band,
                           policy_affine_gap_with_trace_recursion_banded<score_t, trace_t, truncate_score>,
                           policy_affine_gap_with_trace_recursion<score_t, trace_t, truncate_score>>;

    //!\brief Selects the affine gap policy that only computes the score.
    using policy_without_trace_t = std::conditional_t<with_band,
                                                      policy_affine_gap_recursion_banded<score_t, truncate_score>,
                                                      policy_affine_gap_recursion<score_t, truncate_score>>;

public:
    //!\brief The configured recursion policy.
    using type = std::conditional_t<with_trace, policy_with_trace_t, policy_without_trace_t>;
};

/*!\brief Helper type alias for selecting the type of the optimum tracker policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_recursion_policy
 */
template <typename algorithm_traits_t>
using select_recursion_policy_t = typename select_recursion_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
