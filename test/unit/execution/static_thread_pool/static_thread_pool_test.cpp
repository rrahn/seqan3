// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/execution/type_traits/connect_result.hpp>
#include <seqan3/execution/static_thread_pool/static_thread_pool.hpp>
#include <seqan3/execution/sender/sync_wait.hpp>

struct hello_receiver
{

    void set_value()
    {
        std::cout << "hello from thread: " << std::this_thread::get_id() << '\n';
    }

    void set_error(std::exception_ptr eptr)
    {
        try
        {
            if (eptr)
                std::rethrow_exception(eptr);
        }
        catch(const std::exception& e)
        {
            std::cout << "Caught exception \"" << e.what() << "\"\n";
        }
    }

    void set_done()
    {
        // do nothing
    }
};

TEST(static_thread_pool, basic_test)
{
    seqan3::execution::static_thread_pool pool{4};

    auto pool_scheduler = pool.scheduler();
    auto pool_sender = pool_scheduler.schedule();
    auto operation = std::move(pool_sender).connect(hello_receiver{});
    std::cout << "main thread: " << std::this_thread::get_id() << '\n';
    operation.start(); // now we execute on the thread pool.
    seqan3::execution::sync_wait(pool.stop()); // Notifies that the work is done on the pool.
}

// template <typename scheduler_t>
// auto make_task_graph(scheduler_t && scheduler)
// {
//     namespace exec = seqan3::execution;
//     // Now we like to define some kind of task graph that represents a sequence of tasks!
//     // So if the first one is executed the second one will be executed and then the third and so on.
//     // Where lives the operation?
//     auto log_value = [] (int v)
//     {
//         std::cout << "value: " << v "\n";
//         return v;
//     };

//     auto inc_value = [] (int v)
//     {
//         return v + 1;
//     };

//     // So in general we generate a pipeline of nested senders
//     // And only the very last one is connected with a receiver that is supposed to do something with the result.
//     // For this we need to be able to propagate the values generated from one sender to another
//     // Another problem is that the respective sender that gets the values needs to be submitted to the
//     // senders associated context once it has a resolved dependency.


//     // From the pipeline

//     // Like a bulk_execute we want for every task to create a sender to the execution context
//     auto task0 = exec::sequence(exec::on(exec::just(0), scheduler),
//     auto task1 = exec::on(exec::then(exec::transform(std::move(task0), inc_value), log_value)), scheduler);
//     auto task2 = exec::on(exec::then(exec::transform(std::move(task1), inc_value), log_value)), scheduler);
//     auto task3 = exec::on(exec::then(exec::transform(std::move(task2), inc_value), log_value)), scheduler);
//     auto task4 = exec::on(exec::then(exec::transform(std::move(task3), inc_value), log_value)), scheduler);
//     // Now what?
//     // So these define the tasks.
//     // We want to connect each task with an execution context

//     // Simple task graph:
//     // a) create a chain of operations iteratively
//     //   * execute A; A.B.set_value(...)
//     //      * execute B; B.C.set_value(...);
//     //         * execute C; C.print.set_value(...);
//     //             sync_wait()
// }

// TEST(static_thread_pool, task_graph)
// {
//     seqan3::execution::static_thread_pool pool{4};

//     auto task_graph_sender = make_task_graph(pool.scheduler());

//     auto callback = [] (auto value)
//     {
//         std::cout << "final value: " << value << '\n';
//     };
//     seqan3::execution::sync_wait(seqan3::exceution::then(task_graph_sender, callback));
// }
