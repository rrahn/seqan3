// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>
#include <tuple>
#include <seqan3/std/type_traits>

#include <seqan3/execution/cpo/connect.hpp>
#include <seqan3/execution/sender/just.hpp>

template <typename ...values_t>
struct test_just_receiver
{
    std::tuple<std::remove_cvref_t<values_t>...> expected_values;

    test_just_receiver(values_t ...values) : expected_values{std::move(values)...}
    {}

    template <typename ...actual_values_t>
    void set_value(actual_values_t && ...actual_values)
    {
        EXPECT_EQ(std::tie(actual_values...), expected_values);
    }

    template <typename error_t>
    void set_error(error_t const &) noexcept
    {
        FAIL();
    }

    void set_done() noexcept
    {
        FAIL();
    }
};

TEST(just_sender, start)
{
    using namespace std::literals;
    auto op = seqan3::execution::connect(seqan3::execution::just(0, 1.5, "hello"sv),
                                         test_just_receiver{0, 1.5, "hello"sv});
    op.start();
}
