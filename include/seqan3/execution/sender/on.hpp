// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::on.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <optional>

#include <seqan3/execution/cpo/connect.hpp>
#include <seqan3/execution/cpo/schedule.hpp>
#include <seqan3/execution/cpo/start.hpp>

namespace seqan3::execution::detail
{

/*
Constructs a sender s2 such that when connect is called with some receiver output_receiver as
execution::connect(s2, output_receiver) resulting in an operation_state os which is stored as a subobject of the
parent operation_state:
    Constructs a receiver, r and passes r to execution::connect(s, r) resulting in an operation state ros, which is
    stored as a subobject of os such that:
    When set_value, set_error or set_done is called on r, the parameter is copied and stored as a subobject of a receiver r2 and execution::connect(execution::schedule(sch), std::move(r2)) results in an operation_state os2 which is stored as a subobject of os such that:
        When set_value is called on r2, os2’s destructor will be called, the stored value is forwarded to output_receiver on the appropriate choice of set_value, set_error or set_done to match the operation performed on r.
        When set_error or set_done is called on r2 the parameters propagate to output_receiver.
    If connect throws, the resulting exception is forwarded to execution::set_error(output_receiver).
    The destructor of ros is called.
    If connect throws, the resulting exception is forwarded to execution::set_error(output_receiver).
    Calls execution::start(os2).

When execution::start is called on os, call execution::start(r
*/

template <typename value_types, typename error_types>
struct forward_completion_signal
{
    value_types values; // a variant over tuple types
    error_types error;  // a variant over exception types
    bool is_succeeded{false};
    bool is_errored{false};
    bool is_cancelled{false}; // always a bool.
};

template <typename operation_state_t, typename sender_t>
class on_context_receiver
{
    using value_types_t = typename sender_traits<sender_t>::value_types;
    using error_types_t = typename sender_traits<sender_t>::error_traits;

    constexpr sends_done = sender_traits<sender_t>::sends_done;

    forward_completion_signal<value_types_t, error_types_t> fwd_signal;

    operation_state_t & parent_state;

public:
    on_context_receiver() = delete;
    on_context_receiver(on_context_receiver const &) = default;
    on_context_receiver(on_context_receiver &&) = default;
    on_context_receiver & operator=(on_context_receiver const &) = default;
    on_context_receiver & operator=(on_context_receiver &&) = default;
    ~on_context_receiver() = default;

    explicit on_context_receiver(operation_state_t & state,
                                 forward_completion_signal && signal) : parant_state{state},
                                                                        fwd_signal{std::move(signal)}
    {}

    void set_value()
    {
//    When set_value is called on r2, os2’s destructor will be called, the stored value is forwarded to
//    output_receiver on the appropriate choice of set_value, set_error or set_done to match the operation performed on r.
        parent_state.on_context_operation.~on_context_operation_state_t(); // is that safe?
        if (fwd_signal.is_succeeded)
            execution::set_value(std::move(parent_state.output_receiver), std::move(fwd_signal.values));
        else if (fwd_signal.is_errored)
            execution::set_error(std::move(parent_state.output_receiver), std::move(fwd_signal.errors));
        else
            execution::set_done(std::move(parent_state.output_receiver));
    }

    void set_error(std::exception_ptr && error) noexcept
    {
        execution::set_error(std::move(parent_state.output_receiver), std::move(error));
    }

    void set_done() noexcept
    {
        execution::set_done(std::move(parent_state.output_receiver));
    }
};

// Returns a sender that ensures that sender is started on the execution context associated with the specified scheduler.
// The sender is executed with a receiver that customises the get_scheduler query to return the specified scheduler.
// The default implementation schedules the call to connect() and subsequent start() onto an execution context
// associated with scheduler using the schedule(scheduler) operation.
// If schedule(scheduler) completes with set_done() or set_error() then the on() operation completes with that
// signal and never starts executing sender.


// Constructs a sender s2 such that when connect is called with some receiver output_receiver as execution::connect(s2, output_receiver) resulting in an operation_state os which is stored as a subobject of the parent operation_state:
//     Constructs a receiver, r and passes r to execution::connect(s, r) resulting in an operation state ros, which is stored as a subobject of os such that:
//     When set_value, set_error or set_done is called on r, the parameter is copied and stored as a subobject of a receiver r2 and execution::connect(execution::schedule(sch), std::move(r2)) results in an operation_state os2 which is stored as a subobject of os such that:

//         When set_error or set_done is called on r2 the parameters propagate to output_receiver.
//     If connect throws, the resulting exception is forwarded to execution::set_error(output_receiver).
//     The destructor of ros is called.
//     If connect throws, the resulting exception is forwarded to execution::set_error(output_receiver).
//     Calls execution::start(os2).
// When execution::start is called on os, call execution::start(ros).


// on(sender, scheduler) -> on_sender // holds the sender and the scheduler
//     connect(on_sender, output_receiver) -> on_operation_state // if errors, the caller has to handle the error; what happens inside?
//         on_operation_state.start()  // execute operation on the execution context provided by scheduler!
//             // on_operation_state.scheduled_operation = connect(s, out_receiver) && start

template <typename operation_state_t, typename scheduler_t>
class switch_context_receiver
{
private:

    scheduler_t scheduler;
    operation_state_t & state;
public:
    // We schedule this operation here.
    switch_context_receiver(operation_state_t & op_state, scheduler_t && scheduler) :
        state{op_state},
        scheduler{std::move(scheduler)}
    {}

//  When set_value, set_error or set_done is called on r,
//  the parameter is copied and stored as a subobject of a receiver r2 and execution::connect(execution::schedule(sch), std::move(r2))
//  results in an operation_state os2 which is stored as a subobject of os such that:
//    When set_value is called on r2, os2’s destructor will be called, the stored value is forwarded to output_receiver on the appropriate choice of set_value, set_error or set_done to match the operation performed on r.
//    When set_error or set_done is called on r2 the parameters propagate to output_receiver.

    template <typename ...values_t>
    void set_value(values_t && ...values)
    {
        // We forward the values to a second receiver.
        // We need to know the types of the sender
        // This means we need to use sender_traits to figure out what is the type of the receiver.
        on_context_receiver{std::forwards<values_t>(values)...};
        state.
    }

    template <typename error_t>
    void set_error(error_t && error) noexcept
    {
        on_context_receiver{std::forwards<values_t>(values)...};
        execution::set_error(std::move(state.output_receiver), std::current_exception());
    }

    void set_done() noexcept
    {
        execution::set_done(std::move(state.output_receiver));
    }

private:

};

template <typename sender_t, typename scheduler_t>
class on_sender
{
    sender_t sender{};
    scheduler_t scheduler{};

    template <typename receiver_t>
    struct on_operation_state
    {
        receiver_t output_receiver{};

        using sender_operation_state_t = execution::connect_state_t<sender_t,
                                                                    switch_context_receiver<on_operation_state>>;

        using scheduler_sender_t = decltype(execution::schedule(scheduler));
        using on_context_operation_state_t = execution::connect_state_t<scheduler_sender_t,
                                                                        on_context_receiver<on_operation_state,
                                                                                            sender_t>>;

        sender_operation_state_t sender_operation;
        on_context_operation_state_t on_context_operation;

        on_operation_state(sender_t && sender, scheduler_t && scheduler, receiver_t && receiver) :
            scheduler{std::move(scheduler)},
            output_receiver{std::move(receiver)}
        {
            try
            {
                sender_operation = execution::connect(sender, switch_context_receiver{*this, std::move(scheduler)});
            }
            catch (...)
            {
                execution::set_error(std::move(output_receiver), std::current_exception());
            }
        }

        void start() noexcept
        {
            sender_operation.start();
        }
    };

public:

    on_sender() = delete;
    on_sender(on_sender const &) = delete;
    on_sender(on_sender &&) = default;
    on_sender(on_sender const &) = delete;
    on_sender(on_sender &&) = default;
    on_sender(sender_t && sender, scheduler_t scheduler)
        noexcept(std::is_move_constructible_v<sender_t> &&
                 std::is_move_constructible_v<scheduler_t>)
    : sender{std::move(sender)}, scheduler{std::move(scheduler)}
    {}

    template <typename receiver_t>
    on_operation_state connect(receiver_t && output_receiver) &&
    {
        return on_operation_state{std::move(sender), std::move(scheduler), std::forward<receiver_t>(output_receiver)};
    }
};

}  // namespace seqan3::execution::detail

namespace seqan3::execution
{

inline constexpr auto on = [] (auto sender, auto scheduler)
{
    return detail::on_sender{std::move(sender), std::move(scheduler)};
};
}  // namespace seqan3::execution
