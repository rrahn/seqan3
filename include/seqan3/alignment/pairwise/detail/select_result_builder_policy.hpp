// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_result_builder_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_alignment_result_builder.hpp>

namespace seqan3::detail
{

/*!\brief Selects the result builder policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 */
template <typename algorithm_traits_t>
struct select_result_builder_policy
{
    //!\brief The configured alignment result builder policy.
    using type = policy_alignment_result_builder<typename algorithm_traits_t::alignment_result_type>;
};

/*!\brief Helper type alias for selecting the type of the result builder policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_result_builder_policy
 */
template <typename algorithm_traits_t>
using select_result_builder_policy_t = typename select_result_builder_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
