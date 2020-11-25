// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_element_affine_cell.hpp.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <tuple>
#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
template <std::semiregular score_t>
struct alignment_matrix_element_affine_cell
{
    using value_type = score_t;
    using column_value_type = std::pair<value_type, value_type>;
    using vertical_value_type = score_t;

    constexpr static column_value_type initialise_column_value(score_t init) noexcept
    {
        return {init, init};
    }

    constexpr static vertical_value_type initialise_vertical_value(score_t init) noexcept
    {
        return init;
    }
};

}  // namespace seqan3::detail
