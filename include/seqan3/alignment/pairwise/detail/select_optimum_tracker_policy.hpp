// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_optimum_tracker_traits.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker.hpp>
#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker_simd.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?

template <typename algorithm_traits_t, bool is_vectorised>
struct select_optimum_policy
{
private:
    static constexpr bool local_or_no_band = algorithm_traits_t::is_local || !algorithm_traits_t::is_banded;

    using updater_t = std::conditional_t<local_or_no_band, max_score_updater, max_score_banded_updater>;

    using optimum_traits_t = optimum_tracker_traits<typename algorithm_traits_t::score_type,
                                                    typename algorithm_traits_t::matrix_coordinate_type>;

public:
    using type = policy_optimum_tracker<optimum_traits_t, updater_t>;
};

template <typename algorithm_traits_t>
struct select_optimum_policy<algorithm_traits_t, true>
{
    using optimum_tracker_simd_traits_t =
        optimum_tracker_simd_traits<typename algorithm_traits_t::score_type,
                                    typename algorithm_traits_t::original_score_type,
                                    typename algorithm_traits_t::matrix_coordinate_type,
                                    typename algorithm_traits_t::matrix_index_type>;

public:
    using type = policy_optimum_tracker_simd<optimum_tracker_simd_traits_t, max_score_updater_simd_global>;
};

template <typename algorithm_traits_t>
using select_optimum_policy_t = typename select_optimum_policy<algorithm_traits_t,
                                                               algorithm_traits_t::is_vectorised>::type;
}  // namespace seqan3::detail
