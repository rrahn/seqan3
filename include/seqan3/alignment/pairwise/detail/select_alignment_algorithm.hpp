// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/pairwise_alignment_kernel.hpp>

#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm_banded.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{
/*!\brief Selects the pairwise alignment algorithm from the algorithm type traits and the given policies.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 * \tparam policies_t A template parameter pack with the supported policies.
 */
template <typename algorithm_traits_t,
          typename recursion_policy_t,
          typename scoring_scheme_policy_t,
          typename tracker_policy_t,
          typename ...policies_t>
struct select_alignment_algorithm
{
private:
    //!\brief The configured score type.
    using score_t = typename algorithm_traits_t::score_type;
    //!\brief The configured result type.
    using result_t = typename algorithm_traits_t::alignment_result_type;

    using state_t = pairwise_alignment_state<recursion_policy_t, scoring_scheme_policy_t, tracker_policy_t>;

public:
    //!\brief The alignment algorithm.
    using type = lazy_conditional_t<algorithm_traits_t::is_banded,
                                    lazy<pairwise_alignment_algorithm_banded, score_t, result_t, state_t, recursion_policy_t, scoring_scheme_policy_t, tracker_policy_t, policies_t...>,
                                    lazy<pairwise_alignment_algorithm, score_t, result_t, state_t, recursion_policy_t, scoring_scheme_policy_t, tracker_policy_t, policies_t...>>;
};

/*!\brief Helper type alias for selecting the type of the pairwise alignment algorithm.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_alignment_algorithm
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 * \tparam policies_t A template parameter pack with the supported policies.
 */
template <typename algorithm_traits_t, typename ...policies_t>
using select_alignment_algorithm_t = typename select_alignment_algorithm<algorithm_traits_t, policies_t...>::type;

}  // namespace seqan3::detail
