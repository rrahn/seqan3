// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

//!\cond DEV
/*!\file
 * \brief Provides seqan3::align_cfg::simd_gap.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{

/*!\brief Sets the scoring scheme for the alignment algorithm.
 * \ingroup configuration
 * \tparam scoring_scheme_t The type of the scoring scheme. Must satisfy seqan3::ScoringScheme.
 *
 * \details
 *
 * The scoring scheme allows to specify how two symbols of an alphabet are scored inside of the alignment algorithm.
 * The scheme depends on the alphabet type of the passed sequences and must be chosen accordingly.
 * During the configuration of the pairwise alignment algorithm a static assert is triggered if the scoring scheme
 * is not compatible with the given alphabet types. Accordingly, this configuration cannot
 * be defaulted since it depends on the sequences and must be given as a minimal configuration.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/minimal_alignment_config.cpp
 */
template <typename scoring_scheme_t>
struct simd_gap : public pipeable_config_element<simd_gap<scoring_scheme_t>, scoring_scheme_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::simd_gap};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::simd_gap
 * \{
 */

//!\brief Deduces the scoring scheme type from the constructor argument.
template <typename scheme_t>
simd_gap(scheme_t) -> simd_gap<remove_cvref_t<scheme_t>>;
//!\}

} // namespace seqan3::align_cfg
//!\endcond
