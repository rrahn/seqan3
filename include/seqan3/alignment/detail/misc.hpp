// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides miscellaneous utility functions.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/concept.hpp>

namespace seqan3::detail
{

template <typename target_t, typename source_t>
constexpr auto to_simd_if([[maybe_unused]] source_t const t) noexcept
{
    if constexpr (Simd<target_t>)
    {
        return simd::fill<target_t>(t);
    }
    else if constexpr (detail::decays_to_ignore_v<std::remove_reference_t<target_t>>)
    {
        return std::ignore;
    }
    else
    {
        static_assert(std::ConvertibleTo<source_t, target_t>, "The source and target type must be convertible.");
        return static_cast<target_t>(t);
    }
}
} // namespace seqan3::detail
