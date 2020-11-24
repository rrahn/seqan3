// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_element_identity.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace seqan3::detail
{
template <typename trace_t>
struct alignment_matrix_element_identity
{
public:
    using value_type = trace_t;
    using storage_value_type = value_type;
    using element_type = storage_value_type;
    using element_reference = std::add_lvalue_reference_t<storage_value_type>;

    static constexpr bool returns_element_proxy = false;

    constexpr element_reference make_element(storage_value_type & storage_value) const
    {
        return storage_value;
    }

    constexpr storage_value_type initialise(value_type const value) const noexcept
    {
        return value;
    }
};
}  // namespace seqan3::detail
