// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/concepts>
#include <functional>

#include <seqan3/search/configuration/on_result.hpp>

// -----------------------------------------------------------------------------
// Test capturing various callbacks
// -----------------------------------------------------------------------------

TEST(search_cfg_on_result, with_capture_less_lambda)
{
    seqan3::search_cfg::on_result on_result_cfg{[] (auto && result) { return result; }};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}

TEST(search_cfg_on_result, with_capturing_lambda)
{
    int global_result = 0;
    seqan3::search_cfg::on_result on_result_cfg{[&] (auto && result) { global_result = result; }};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ(global_result, 0);
    std::invoke(on_result_cfg.callback, 10);
    EXPECT_EQ(global_result, 10);
}

int my_free_function(int v)
{
    return v;
}

TEST(search_cfg_on_result, with_free_function)
{
    seqan3::search_cfg::on_result on_result_cfg{my_free_function};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}

struct my_function_object
{
    template <typename t>
    t operator()(t && v)
    {
        return std::forward<t>(v);
    }
};

TEST(search_cfg_on_result, with_function_object)
{
    seqan3::search_cfg::on_result on_result_cfg{my_function_object{}};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}
