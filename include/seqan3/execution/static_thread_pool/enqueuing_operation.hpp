// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool_enqueuing_operation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/execution/static_thread_pool/operation_base.hpp>

namespace seqan3::execution
{

template <typename queue_t, typename receiver_t>
class static_thread_pool_enqueuing_operation : public static_thread_pool_operation_base
{
private:
    //!\brief The structure to enqueue a new operation.
    queue_t & queue;
    //!\brief The receiver to invoke.
    receiver_t receiver;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    static_thread_pool_enqueuing_operation() = default; //!< Default.
    static_thread_pool_enqueuing_operation(static_thread_pool_enqueuing_operation const &) = delete; //!< Delete.
    static_thread_pool_enqueuing_operation(static_thread_pool_enqueuing_operation &&) = delete; //!< Delete.
    static_thread_pool_enqueuing_operation & operator=(static_thread_pool_enqueuing_operation const &) = delete; //!< Delete.
    static_thread_pool_enqueuing_operation & operator=(static_thread_pool_enqueuing_operation &&) = delete; //!< Delete.
    ~static_thread_pool_enqueuing_operation() = default; //!< Default.

    static_thread_pool_enqueuing_operation(queue_t & queue, receiver_t receiver) :
        queue{queue},
        receiver{std::move(receiver)}
    {}
    //!\}

    //!\brief Calls the receiver's set value function.
    void set_value() noexcept override
    {
        // Which of the two interfaces need to be called?
        // How am I providing the necessary data to be moved to the next execution?
        try
        {
            receiver.set_value();
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
    void start() noexcept
    {
        [[maybe_unused]] contrib::queue_op_status push_status = queue.wait_push(this);
        assert(push_status == contrib::queue_op_status::success);
    }
};

}  // namespace seqan3::execution
