// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides execution seqan3::execution::connect.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::execution::cpo
{

// poison pill
void connect() = delete;

class connect_fn
{
public:
    template <typename sender_t, typename receiver_t>
    auto operator()(sender_t && sender, receiver_t && receiver) const noexcept
    {
        return std::move(sender).connect(std::forward<receiver_t>(receiver));
    }
};
} // seqan3::execution::cpo

namespace seqan3::execution
{

inline constexpr cpo::connect_fn connect{};

}  // namespace seqan3::execution
