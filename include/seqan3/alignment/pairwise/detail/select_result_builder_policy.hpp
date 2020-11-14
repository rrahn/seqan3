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
// TODO: Can we also work with forward declarations here?

template <typename algorithm_traits_t>
struct select_result_builder_policy
{
private:
    using result_builder_traits_t =
        alignment_result_builder_traits<typename algorithm_traits_t::alignment_result_type>;

public:
    //!\brief The configured recursion policy.
    using type = policy_alignment_result_builder<result_builder_traits_t>;
};

template <typename algorithm_traits_t>
using select_result_builder_policy_t = typename select_result_builder_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
