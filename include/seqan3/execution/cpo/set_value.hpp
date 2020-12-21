// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides execution seqan3::execution::set_value.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::execution::cpo
{

// poison pill
void set_value() = delete;

class set_value_fn
{
public:
    template <typename receiver_t, typename ...values_t>
    void operator()(receiver_t && receiver, values_t && ...values) const noexcept
    {
        std::forward<receiver_t>(receiver).set_value(std::forward<values_t>(values)...);
    }
};
} // seqan3::execution::cpo

namespace seqan3::execution
{

inline constexpr cpo::set_value_fn set_value{};

}  // namespace seqan3::execution
