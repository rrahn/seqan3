// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool.hpp.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <thread>
#include <vector>

#include <seqan3/contrib/parallel/buffer_queue.hpp>

#include <seqan3/execution/static_thread_pool/operation_base.hpp>
#include <seqan3/execution/static_thread_pool/enqueuing_sender.hpp>
#include <seqan3/execution/static_thread_pool/terminating_sender.hpp>
#include <seqan3/utility/parallel/detail/spin_delay.hpp>

namespace seqan3::execution
{

template <typename, typename>
class static_thread_pool_terminating_operation;

class static_thread_pool
{
private:

    template <typename, typename>
    friend class static_thread_pool_terminating_operation;

    //!\brief The operation queue type.
    using operation_queue_t = contrib::dynamic_buffer_queue<static_thread_pool_operation_base *>;

    //!\brief The thread pool enabling concurrency.
    std::vector<std::thread> worker_pool{};
    operation_queue_t operation_queue{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    static_thread_pool() = delete;
    static_thread_pool(static_thread_pool const &) = delete;
    static_thread_pool(static_thread_pool &&) = delete;
    static_thread_pool & operator=(static_thread_pool const &) = delete;
    static_thread_pool & operator=(static_thread_pool &&) = delete;

    static_thread_pool(size_t const thread_count)
    {
        for (size_t i = 0; i < thread_count; ++i)
        {
            worker_pool.emplace_back([this] ()
            {
                while (true)
                {
                    static_thread_pool_operation_base * operation = nullptr;
                    if (contrib::queue_op_status pop_status = operation_queue.wait_pop(operation);
                        pop_status == contrib::queue_op_status::closed)
                    {
                        break; // no more work is comming so stop spinning and release the thread resource.
                    }
                    else if (pop_status != contrib::queue_op_status::success)
                    {
                        assert(pop_status != contrib::queue_op_status::closed);
                        break;
                        // todo: better error handling
                        // signal interrupt to all remaining operations
                        // operation_queue.close(); -> no more work can come
                        // producer thread needs to be aware of this
                        // TODO: think about reopen the queue, so we should not free thread resource.
                    }
                    // The set value is valid.
                    // assert(pop_status == contrib::queue_op_status::success);
                    assert(operation != nullptr); // TODO: better error handling?

                    operation->set_value(); // notify the actual operation to execute and send to receiver.
                }
            });
        }
    }

    ~static_thread_pool()
    {
        close_and_wait();
        // Make sure all threads are joined.
        for (std::thread & worker : worker_pool)
            if (worker.joinable())
                worker.join();
    }
    //!\}

    struct static_thread_pool_scheduler
    {
        operation_queue_t & operation_queue;

        auto schedule() noexcept
        {
            return static_thread_pool_enqueuing_sender{operation_queue};
        }
    };

    static_thread_pool_scheduler scheduler() noexcept
    {
        return static_thread_pool_scheduler{operation_queue};
    }

    auto stop() noexcept
    {
        return static_thread_pool_terminating_sender{*this};
    }
private:

    void close_and_wait() noexcept
    {
        operation_queue.close();

        detail::spin_delay delay{};
        while (!operation_queue.is_empty())
            delay.wait();

        assert(operation_queue.is_empty());
    }
};
}  // namespace seqan3::execution
