// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides execution seqan3::execution::schedule.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::execution::cpo
{

// poison pill
void schedule() = delete;

class schedule_fn
{
public:
    template <typename scheduler_t>
    auto operator()(scheduler_t && scheduler) const noexcept
    {
        return std::forward<scheduler_t>(scheduler).schedule();
    }
};
} // seqan3::execution::cpo

namespace seqan3::execution
{

inline constexpr cpo::schedule_fn schedule{};

}  // namespace seqan3::execution
