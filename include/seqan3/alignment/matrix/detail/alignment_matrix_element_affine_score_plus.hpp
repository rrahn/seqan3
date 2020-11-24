// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_element_affine_score_plus.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>
#include <tuple>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
template <typename score_t>
struct alignment_matrix_element_affine_score_plus
{
private:
    std::pair<score_t, score_t> vertical_score{};
public:

    using value_type = score_t;
    using storage_value_type = std::pair<value_type, value_type>;
    using element_type = affine_cell<std::tuple<value_type, value_type, value_type, value_type>>;
    using element_reference = affine_cell<std::tuple<value_type &, value_type &, value_type &, value_type &>>;

    static constexpr bool returns_element_proxy = true;

    constexpr element_reference make_element(storage_value_type & storage_value)
    {
        return {storage_value.first, storage_value.second, vertical_score.first, vertical_score.second};
    }

    constexpr storage_value_type initialise(value_type const value) noexcept
    {
        vertical_score = {value, value};
        return {value, value};
    }
};
}  // namespace seqan3::detail
