// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::wait.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <condition_variable>
#include <mutex>

namespace seqan3::execution::detail
{

template <typename synchronised_state_t>
class sync_wait_receiver
{

    synchronised_state_t & sync_state;
public:

    explicit sync_wait_receiver(synchronised_state_t & state) noexcept : sync_state{state}
    {}

    sync_wait_receiver() = delete;
    sync_wait_receiver(sync_wait_receiver const &) = delete;
    sync_wait_receiver(sync_wait_receiver &&) noexcept = default;
    sync_wait_receiver & operator=(sync_wait_receiver const &) = delete;
    sync_wait_receiver & operator=(sync_wait_receiver &&) noexcept = default;
    ~sync_wait_receiver() = default;

    // template <typename value_t>
    // void set_value(value_t && value) noexcept
    // {
    //     {
    //         std::scoped_lock notify_lock{sync_state.condition_mutex};
    //         sync_state.is_done = true;
    //     }
    //     sync_state.condition.notify_one();
    // }

    void set_value() noexcept
    {
        notify_sender();
    }

    void set_error(std::exception_ptr error) noexcept
    {
        sync_state.error = error;
        notify_sender();
    }

    void set_done() noexcept
    {
        sync_state.is_cancelled = true;
        notify_sender();
    }
private:

    void notify_sender() noexcept
    {
        {
            std::scoped_lock notify_lock{sync_state.condition_mutex};
            sync_state.is_done = true;
        }
        sync_state.condition.notify_one();
    }
};

} // namespace seqan3::execution::detail

namespace seqan3::execution
{

inline constexpr auto sync_wait = [] (auto && sender) -> void
{
    struct synchronised_state_t
    {
        std::mutex condition_mutex;
        std::condition_variable condition;
        std::optional<std::exception_ptr> error;
        bool is_done{false};
        bool is_cancelled{false};
    } sync_state;

    auto operation = std::move(sender).connect(detail::sync_wait_receiver<synchronised_state_t>{sync_state});
    operation.start(); // submit work.

    // Waiting until work is done.
    std::unique_lock wait_lock(sync_state.condition_mutex);
    sync_state.condition.wait(wait_lock, [&] { return sync_state.is_done; });

    // Check synchronised state for anomalies.
    if (sync_state.error.has_value())
        std::rethrow_exception(*sync_state.error);

    if (sync_state.is_cancelled)
        std::terminate();
};

}  // namespace seqan3::execution
