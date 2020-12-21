// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

template <typename sender1_t, typename sender2_t>
struct notify_receiver
{


    template <typename ...values_t>
    void set_done(values_t && ...values)
    {
        // signalling to the next senders that they can be started.
        // we simply call submit and let the thread context handle the memory?

    }

    template <typename error_t>
    void set_error(error_t && error) noexcept
    {
        // Now we have to handle the error here.
    }

    template <typename error_t>
    void set_done() noexcept
    {
        // now, what to do?
        // We should clean any open resources
        // and call done on the original thing -> that is no more work will be executed.
        // Which leaves us with the question, what happens with work that is still open yet?
        //
    }
};

TEST(task_graph, create)
{
    // We need dynamic storage of the sender
    struct sender_t
    {
        schedule(scheduler) | then(compute_algorithm);
    } sender;

    auto sender = scheduler.schedule();

    auto algo = [](int v) { return v + 1; };
    auto compute_block = sender | execution::then(algo); // So every block is increasing a value:

    // First step create graph:
    std::vector<decltype(compute_block)> graph{};
    for (unsigned i = 0; i < 10; ++i)
        graph.emplace_back(execution::then(scheduler.schedule(), algo));

    // Second step connect graph
    for (int i = 1; i < 10; ++i)
        (graph[i-1], graph[i]) -> so when i-1 is done it will invoke sender i
        // where is this stored?



    auto result = execution::sync_wait(execution::as_graph(graph.front(), graph.back()));

    // as_graph returns a sender that starts the execution on the first node and computes all elements until
    // the last node has been computed, which will call the receiver it connects to with the result.
    // How are we modeling the operation states in between?
    // How doe we know to connect a new node?
    // Where is the dependency comming from?
    // That must be set before?
}
