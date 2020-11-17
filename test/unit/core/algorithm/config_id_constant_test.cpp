// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "configuration_mock.hpp"

#include <seqan3/core/algorithm/config_id_constant.hpp>

TEST(configuration_utility, config_id_constant)
{
    seqan3::detail::config_id_constant<test_algo_id> config_ids{test_algo_id::bar_id,
                                                                test_algo_id::bax_id,
                                                                test_algo_id::foobar_id};

    EXPECT_TRUE(config_ids.contains(test_algo_id::bar_id));
    EXPECT_TRUE(config_ids.contains(test_algo_id::bax_id));
    EXPECT_TRUE(config_ids.contains(test_algo_id::foobar_id));
    EXPECT_FALSE(config_ids.contains(test_algo_id::foo_id));
}

TEST(configuration_utility, config_id_constant_constexpr)
{
    constexpr seqan3::detail::config_id_constant<test_algo_id> config_ids{test_algo_id::bar_id,
                                                                          test_algo_id::bax_id,
                                                                          test_algo_id::foobar_id};

    constexpr bool contains_bar_id = config_ids.contains(test_algo_id::bar_id);
    constexpr bool contains_bax_id = config_ids.contains(test_algo_id::bax_id);
    constexpr bool contains_foobar_id = config_ids.contains(test_algo_id::foobar_id);
    constexpr bool contains_foo_id = config_ids.contains(test_algo_id::foo_id);

    EXPECT_TRUE(contains_bar_id);
    EXPECT_TRUE(contains_bax_id);
    EXPECT_TRUE(contains_foobar_id);
    EXPECT_FALSE(contains_foo_id);
}
