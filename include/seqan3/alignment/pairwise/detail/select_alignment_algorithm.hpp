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

#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm_banded.hpp>
#include <seqan3/core/type_traits/lazy.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?
template <typename algorithm_traits_t, typename ...policies_t>
struct select_alignment_algorithm
{
private:
    //!\brief Selects either the banded or the unbanded alignment algorithm based on the given traits type.
    using algo_traits_t = algorithm_traits<typename algorithm_traits_t::score_type,
                                           typename algorithm_traits_t::original_score_type,
                                           typename algorithm_traits_t::alignment_result_type>;

public:
    //!\brief The configured recursion policy.
    using type = lazy_conditional_t<algorithm_traits_t::is_banded,
                                    lazy<pairwise_alignment_algorithm_banded, algo_traits_t, policies_t...>,
                                    lazy<pairwise_alignment_algorithm, algo_traits_t, policies_t...>>;
};

template <typename algorithm_traits_t, typename ...policies_t>
using select_alignment_algorithm_t = typename select_alignment_algorithm<algorithm_traits_t, policies_t...>::type;

}  // namespace seqan3::detail
