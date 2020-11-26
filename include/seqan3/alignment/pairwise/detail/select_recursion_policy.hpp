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
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion_opt.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion_banded.hpp>
#include <seqan3/core/type_traits/lazy.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?

template <typename algorithm_traits_t,
          template <typename ...> typename policy_score,
          template <typename ...> typename policy_with_trace>
struct select_recursion_policy
{
private:
    static constexpr bool with_trace = algorithm_traits_t::requires_trace_information;

    using recursion_policy_traits_t =
        affine_gap_recursion_traits<typename algorithm_traits_t::score_type,
                                    typename algorithm_traits_t::original_score_type,
                                    typename algorithm_traits_t::trace_type,
                                    algorithm_traits_t::is_local>;

public:
    //!\brief The configured recursion policy.
    using type = lazy_conditional_t<with_trace,
                                    lazy<policy_with_trace, recursion_policy_traits_t>,
                                    lazy<policy_score, recursion_policy_traits_t>>;
};

template <typename algorithm_traits_t>
using select_recursion_policy_t =
    // std::conditional_t<algorithm_traits_t::is_banded,
    //                    typename select_recursion_policy<algorithm_traits_t,
    //                                                     policy_affine_gap_recursion_banded,
    //                                                     policy_affine_gap_with_trace_recursion_banded>::type,
                       typename select_recursion_policy<algorithm_traits_t,
                                                        policy_affine_gap_recursion,
                                                        policy_affine_gap_with_trace_recursion_opt>::type;//>;
}  // namespace seqan3::detail
