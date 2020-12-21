// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool_enqueuing_sender.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/execution/static_thread_pool/enqueuing_operation.hpp>

namespace seqan3::execution
{

template <typename operation_queue_t>
class static_thread_pool_enqueuing_sender
{
private:
    // We need to be able to queue here something?
    operation_queue_t & queue;
public:
    explicit static_thread_pool_enqueuing_sender(operation_queue_t & queue) : queue{queue}
    {}
    static_thread_pool_enqueuing_sender(static_thread_pool_enqueuing_sender const &) = delete;
    static_thread_pool_enqueuing_sender(static_thread_pool_enqueuing_sender &&) = default;
    ~static_thread_pool_enqueuing_sender() = default;

    template <typename receiver_t>
    auto connect(receiver_t && receiver) && -> static_thread_pool_enqueuing_operation<operation_queue_t, receiver_t>
    {
        return static_thread_pool_enqueuing_operation{queue, std::forward<receiver_t>(receiver)};
    }
};

}  // namespace seqan3::execution
