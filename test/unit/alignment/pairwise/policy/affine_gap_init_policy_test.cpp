// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/test/simd_utility.hpp>
#include <seqan3/test/pretty_printing.hpp>

class affine_gap_init_policy_mock :
    public seqan3::detail::affine_gap_init_policy<affine_gap_init_policy_mock>
{
public:
    using base_t = seqan3::detail::affine_gap_init_policy<affine_gap_init_policy_mock>;

    using base_t::base_t;
    using base_t::init_origin_cell;
    using base_t::init_column_cell;
    using base_t::init_row_cell;
};

using namespace seqan3;

template <typename t>
struct affine_gap_init_policy_test : public ::testing::Test
{

    t convert(int32_t const v)
    {
        if constexpr (Simd<t>)
            return fill<t>(v);
        else
            return v;
    }
};

using test_types = ::testing::Types<int32_t, simd_type_t<int32_t>>;

TYPED_TEST_CASE(affine_gap_init_policy_test, test_types);

TEST(affine_gap_init_policy, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<affine_gap_init_policy_mock>);
    EXPECT_TRUE(std::is_copy_constructible_v<affine_gap_init_policy_mock>);
    EXPECT_TRUE(std::is_move_constructible_v<affine_gap_init_policy_mock>);
    EXPECT_TRUE(std::is_copy_assignable_v<affine_gap_init_policy_mock>);
    EXPECT_TRUE(std::is_move_assignable_v<affine_gap_init_policy_mock>);
    EXPECT_TRUE(std::is_destructible_v<affine_gap_init_policy_mock>);
}

// How can we make the test also work for simd types?
TYPED_TEST(affine_gap_init_policy_test, init_origin_cell)
{
    std::tuple cell{this->convert(0), this->convert(0), std::ignore};
    std::tuple cache{std::tuple{this->convert(0), this->convert(0), std::ignore}, this->convert(-10), this->convert(-1)};

    affine_gap_init_policy_mock mock{};

    mock.init_origin_cell(std::make_tuple(std::ref(cell), seqan3::alignment_coordinate{}, std::ignore), cache);

    auto & [first, second, trace] = cell;
    (void) trace;
    SIMD_OR_SCALAR_EQ(first, this->convert(0));
    SIMD_OR_SCALAR_EQ(second, this->convert(-10));
    std::tie(first, second, std::ignore) = std::get<0>(cache);
    SIMD_OR_SCALAR_EQ(first, this->convert(0));
    SIMD_OR_SCALAR_EQ(second, this->convert(-10));
    SIMD_OR_SCALAR_EQ(std::get<1>(cache), this->convert(-10));
    SIMD_OR_SCALAR_EQ(std::get<2>(cache), this->convert(-1));
}

TYPED_TEST(affine_gap_init_policy_test, init_column_cell)
{
    std::tuple cell{this->convert(0), this->convert(-10), std::ignore};
    std::tuple cache{std::tuple{this->convert(0), this->convert(-10), std::ignore}, this->convert(-10), this->convert(-1)};

    affine_gap_init_policy_mock mock{};

    mock.init_column_cell(std::make_tuple(std::ref(cell), seqan3::alignment_coordinate{}, std::ignore), cache);

    auto & [first, second, trace] = cell;
    (void) trace;
    SIMD_OR_SCALAR_EQ(first, this->convert(-10));
    SIMD_OR_SCALAR_EQ(second, this->convert(-20));
    std::tie(first, second, std::ignore) = std::get<0>(cache);
    SIMD_OR_SCALAR_EQ(first, this->convert(0));
    SIMD_OR_SCALAR_EQ(second, this->convert(-11));
    SIMD_OR_SCALAR_EQ(std::get<1>(cache), this->convert(-10));
    SIMD_OR_SCALAR_EQ(std::get<2>(cache), this->convert(-1));
}

TYPED_TEST(affine_gap_init_policy_test, init_row_cell)
{
    std::tuple cell{this->convert(0), this->convert(-10), std::ignore};
    std::tuple cache{std::tuple{this->convert(0), this->convert(0), std::ignore}, this->convert(-10), this->convert(-1)};

    affine_gap_init_policy_mock mock{};

    mock.init_row_cell(std::make_tuple(std::ref(cell), seqan3::alignment_coordinate{}, std::ignore), cache);

    auto & [first, second, trace] = cell;
    (void) trace;
    SIMD_OR_SCALAR_EQ(first, this->convert(-10));
    SIMD_OR_SCALAR_EQ(second, this->convert(-11));
    std::tie(first, second, std::ignore) = std::get<0>(cache);
    SIMD_OR_SCALAR_EQ(first, this->convert(0));
    SIMD_OR_SCALAR_EQ(second, this->convert(-20));
    SIMD_OR_SCALAR_EQ(std::get<1>(cache), this->convert(-10));
    SIMD_OR_SCALAR_EQ(std::get<2>(cache), this->convert(-1));
}
