// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool_terminating_sender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/execution/static_thread_pool/terminating_operation.hpp>

namespace seqan3::execution
{

template <typename static_thread_pool_t>
class static_thread_pool_terminating_sender
{
private:

    friend static_thread_pool_t;

    static_thread_pool_t & pool;

    explicit static_thread_pool_terminating_sender(static_thread_pool_t & pool) : pool{pool}
    {}
public:

    static_thread_pool_terminating_sender(static_thread_pool_terminating_sender const &) = delete;
    static_thread_pool_terminating_sender(static_thread_pool_terminating_sender &&) = default;
    ~static_thread_pool_terminating_sender() = default;

    template <typename receiver_t>
    auto connect(receiver_t && receiver) && -> static_thread_pool_terminating_operation<static_thread_pool_t, receiver_t>
    {
        return static_thread_pool_terminating_operation{pool, std::forward<receiver_t>(receiver)};
    }
};

}  // namespace seqan3::execution
