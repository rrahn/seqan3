// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool_terminating_operation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <chrono>

#include <seqan3/execution/static_thread_pool/operation_base.hpp>

namespace seqan3::execution
{

template <typename static_thread_pool_t, typename receiver_t>
class static_thread_pool_terminating_operation : public static_thread_pool_operation_base
{
private:
    //!\brief The structure to enqueue a new operation.
    static_thread_pool_t & pool;
    //!\brief The receiver to invoke.
    receiver_t receiver;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    static_thread_pool_terminating_operation() = default; //!< Default.
    static_thread_pool_terminating_operation(static_thread_pool_terminating_operation const &) = delete; //!< Delete.
    static_thread_pool_terminating_operation(static_thread_pool_terminating_operation &&) = delete; //!< Delete.
    static_thread_pool_terminating_operation & operator=(static_thread_pool_terminating_operation const &) = delete; //!< Delete.
    static_thread_pool_terminating_operation & operator=(static_thread_pool_terminating_operation &&) = delete; //!< Delete.
    ~static_thread_pool_terminating_operation() = default; //!< Default.

    static_thread_pool_terminating_operation(static_thread_pool_t & pool, receiver_t receiver) :
        pool{pool},
        receiver{std::move(receiver)}
    {}
    //!\}

    //!\brief Calls the receiver's set value function.
    void set_value() noexcept override
    {
        using namespace std::chrono_literals;
        // Which of the two interfaces need to be called?
        // How am I providing the necessary data to be moved to the next execution?
        try
        {
            std::this_thread::sleep_for(2s);
            pool.close_and_wait();

            receiver.set_value(); // notify that everything is finished.
        }
        catch (...)
        {
            receiver.set_error(std::current_exception());
        }
    }

    //!\brief Calls the receiver's set value function.
    void set_error() noexcept override
    {
        receiver.set_error(std::current_exception());
    }

    void set_done() noexcept override
    {
        receiver.set_done();
    }

    //!\brief Asynchronously pushes the value.
    void start()
    {
        [[maybe_unused]] contrib::queue_op_status push_status = pool.operation_queue.wait_push(this);
        assert(push_status == contrib::queue_op_status::success);
    }
};

}  // namespace seqan3::execution
