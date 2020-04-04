// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_policy_scoring.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief The CRTP-policy that stores the scoring scheme used for this alignment algorithm.
 * \ingroup alignment_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam scoring_scheme_t The type of the scoring scheme.
 *
 * \remarks The template parameters of this CRTP-policy are selected in the
 *          seqan3::detail::alignment_configurator::configure_scoring_scheme when selecting the alignment for the given
 *          configuration.
 */
template <typename alignment_configuration_t, typename scoring_scheme_t>
class pairwise_alignment_policy_scoring
{
protected:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr pairwise_alignment_policy_scoring() noexcept = default;
    //!\brief Defaulted.
    constexpr pairwise_alignment_policy_scoring(pairwise_alignment_policy_scoring const &) noexcept = default;
    //!\brief Defaulted.
    constexpr pairwise_alignment_policy_scoring(pairwise_alignment_policy_scoring &&) noexcept = default;
    //!\brief Defaulted.
    constexpr pairwise_alignment_policy_scoring & operator=(pairwise_alignment_policy_scoring const &) noexcept = default;
    //!\brief Defaulted.
    constexpr pairwise_alignment_policy_scoring & operator=(pairwise_alignment_policy_scoring &&) noexcept = default;
    //!\brief Defaulted.
    ~pairwise_alignment_policy_scoring() noexcept = default;

    pairwise_alignment_policy_scoring(alignment_configuration_t const & config)
        : scoring_scheme{seqan3::get<align_cfg::scoring>(config).value}
    {}
    //!\}

    //!\brief The scoring scheme used for this alignment algorithm.
    scoring_scheme_t scoring_scheme{};
};

} // namespace seqan3::detail
