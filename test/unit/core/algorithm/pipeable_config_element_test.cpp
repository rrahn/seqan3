// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "configuration_mock.hpp"

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
enum struct incompatible_id : uint8_t
{
    incompatible
};
} // namespace seqan3::detail

class incompatible_config : private seqan3::pipeable_config_element
{
public:

    incompatible_config() = default;
    incompatible_config(incompatible_config const &) = default;
    incompatible_config(incompatible_config &&) = default;
    incompatible_config & operator=(incompatible_config const &) = default;
    incompatible_config & operator=(incompatible_config &&) = default;
    ~incompatible_config() = default;

    static constexpr seqan3::detail::incompatible_id id = seqan3::detail::incompatible_id::incompatible;
};

TEST(pipeable_config_element, pipeable_concepts)
{
    EXPECT_TRUE(seqan3::config_element_specialisation<bar>);
    EXPECT_TRUE(seqan3::config_element_specialisation<bax>);
    EXPECT_TRUE(seqan3::config_element_specialisation<incompatible_config>);
    EXPECT_TRUE((seqan3::config_element_pipeable_with<bar, bax>));
    EXPECT_FALSE((seqan3::config_element_pipeable_with<bar, incompatible_config>));
    EXPECT_FALSE((seqan3::config_element_pipeable_with<incompatible_config, bax>));
}

TEST(pipeable_config_element, two_elements)
{
    bar b1{};
    bax b2{};
    // lvalue | lvalue
    {
        auto cfg = b1 | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | lvalue
    {
        auto cfg = bar{} | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // lvalue | rvalue
    {
        auto cfg = b1 | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = bar{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }
}

TEST(pipeable_config_element, configuration_with_element)
{
    seqan3::configuration<bar> tmp{};
    bax b2{};

    // lvalue | lvalue
    {
        auto cfg = tmp | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | lvalue
    {
        auto cfg = seqan3::configuration<bar>{} | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // lvalue | rvalue
    {
        auto cfg = tmp | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = seqan3::configuration<bar>{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }
}

TEST(pipeable_config_element, element_with_configuration)
{
    seqan3::configuration<bar> tmp{};
    bax b2{};

    // lvalue | lvalue
    {
        auto cfg = b2 | tmp;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax, bar>>));
    }

    // rvalue | lvalue
    {
        auto cfg = bax{} | tmp;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax, bar>>));
    }

    // lvalue | rvalue
    {
        auto cfg = b2 | seqan3::configuration<bar>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax, bar>>));
    }

    // rvalue | rvalue
    {
        auto cfg = bax{} | seqan3::configuration<bar>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax, bar>>));
    }
}

TEST(pipeable_config_element, configuration_with_configuration)
{
    seqan3::configuration<bar> tmp{};
    seqan3::configuration<bax> b2{};
    // lvalue | lvalue
    {
        auto cfg = tmp | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | lvalue
    {
        auto cfg = seqan3::configuration<bar>{} | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // lvalue | rvalue
    {
        auto cfg = tmp | seqan3::configuration<bax>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = seqan3::configuration<bar>{} | seqan3::configuration<bax>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax>>));
    }
}

TEST(pipeable_config_element, special_cases)
{
    // with empty configuration on left hand side
    {
        auto cfg = seqan3::configuration<>{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax>>));
    }

    // with empty configuration on left hand side
    {
        auto cfg = seqan3::configuration<>{} | seqan3::configuration{bax{}};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax>>));
    }

    // with empty configuration on right hand side
    {
        auto cfg = bax{} | seqan3::configuration<>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax>>));
    }

    // with empty configuration on right hand side
    {
        auto cfg = seqan3::configuration{bax{}} | seqan3::configuration<>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bax>>));
    }

    // with two empty configurations
    {
        auto cfg = seqan3::configuration<>{} | seqan3::configuration<>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<>>));
    }
}

TEST(pipeable_config_element, multiple_elements)
{
    seqan3::configuration<bar> tmp{};
    bax b2{};
    {
        auto cfg = foo{} | bar{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        auto cfg = seqan3::configuration<bar>{} | b2 | foo{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<bar, bax, foo>>));
    }
}

TEST(pipeable_config_element, const_config)
{
    seqan3::configuration<foobar<>> const tmp{};

    {
        auto cfg = tmp | foo{} | bar{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<foobar<>, foo, bar>>));
    }
}
