// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides execution seqan3::execution::start.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::execution::cpo
{

// poison pill
void start() = delete;

class start_fn
{
public:
    template <typename operation_state_t>
    auto operator()(operation_state_t && operation_state) const noexcept
    {
        return operation_state.start();
    }
};
} // seqan3::execution::cpo

namespace seqan3::execution
{

inline constexpr cpo::start_fn start{};

}  // namespace seqan3::execution
