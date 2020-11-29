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

namespace seqan3::detail
{
/*!\brief Selects the optimum tracker policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 *
 * \details
 *
 * Configures the optimal tracker policy type based on the configured alignment algorithm traits.
 * The selector differentiates between vectorised and scalar alignment computation. In the scalar mode, the
 * score update function object is configured based on the band and method information of the given traits.
 */
template <typename algorithm_traits_t>
struct select_optimum_policy
{
private:
    //!\brief A flag indicating whether the alignment uses the local method or does not use a band.
    static constexpr bool local_or_no_band = algorithm_traits_t::is_local || !algorithm_traits_t::is_banded;
    //!\brief The type of the update function object depending on the selected alignment method in scalar mode.
    using updater_t = std::conditional_t<local_or_no_band, max_score_updater, max_score_banded_updater>;
    //!\brief The configured score type.
    using score_t = typename algorithm_traits_t::score_type;
    //!\brief The configured matrix coordinate type.
    using coordinate_t = typename algorithm_traits_t::matrix_coordinate_type;

public:
    //!\brief The type of the optimum tracker policy.
    using type =
        lazy_conditional_t<algorithm_traits_t::is_vectorised,
                           lazy<policy_optimum_tracker_simd, score_t, coordinate_t, max_score_updater_simd_global>,
                           lazy<policy_optimum_tracker, score_t, coordinate_t, updater_t>>;
};

/*!\brief Helper type alias for selecting the type of the optimum tracker policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_optimum_policy
 */
template <typename algorithm_traits_t>
using select_optimum_policy_t = typename select_optimum_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
