// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/core/simd/all.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

template <typename t>
struct alignment_optimum_test : public ::testing::Test
{};

using testing_types = ::testing::Types<int32_t, simd_type_t<int32_t>>;
TYPED_TEST_CASE(alignment_optimum_test, testing_types);

TYPED_TEST(alignment_optimum_test, construction)
{
    EXPECT_TRUE((std::is_default_constructible<detail::alignment_optimum<TypeParam>>::value));
    EXPECT_TRUE((std::is_copy_constructible<detail::alignment_optimum<TypeParam>>::value));
    EXPECT_TRUE((std::is_copy_assignable<detail::alignment_optimum<TypeParam>>::value));
    EXPECT_TRUE((std::is_move_constructible<detail::alignment_optimum<TypeParam>>::value));
    EXPECT_TRUE((std::is_move_assignable<detail::alignment_optimum<TypeParam>>::value));
    EXPECT_TRUE((std::is_destructible<detail::alignment_optimum<TypeParam>>::value));
}

TYPED_TEST(alignment_optimum_test, type_deduction)
{
    detail::alignment_optimum def_ao{};
    EXPECT_TRUE((std::is_same_v<decltype(def_ao), detail::alignment_optimum<int32_t>>));

    alignment_coordinate<TypeParam> coordinate{};
    coordinate.first = detail::to_simd_if<TypeParam>(2u);
    coordinate.second = detail::to_simd_if<TypeParam>(3u);
    detail::alignment_optimum ao{detail::to_simd_if<TypeParam>(10), coordinate};
    EXPECT_TRUE((std::is_same_v<decltype(ao), detail::alignment_optimum<TypeParam>>));
}

TYPED_TEST(alignment_optimum_test, access)
{
    alignment_coordinate<TypeParam> coordinate{};
    coordinate.first = detail::to_simd_if<TypeParam>(2u);
    coordinate.second = detail::to_simd_if<TypeParam>(3u);
    detail::alignment_optimum ao{detail::to_simd_if<TypeParam>(10), coordinate};
    SIMD_OR_SCALAR_EQ(ao.score, detail::to_simd_if<TypeParam>(10));
    SIMD_OR_SCALAR_EQ(ao.coordinate.first, detail::to_simd_if<TypeParam>(2));
    SIMD_OR_SCALAR_EQ(ao.coordinate.second, detail::to_simd_if<TypeParam>(3));
}

TYPED_TEST(alignment_optimum_test, max)
{
    alignment_coordinate<TypeParam> coordinate_l{};
    coordinate_l.first = detail::to_simd_if<TypeParam>(1u);
    coordinate_l.second = detail::to_simd_if<TypeParam>(4u);
    detail::alignment_optimum aol{detail::to_simd_if<TypeParam>(9), coordinate_l};
    SIMD_OR_SCALAR_EQ(aol.score, detail::to_simd_if<TypeParam>(9));
    SIMD_OR_SCALAR_EQ(aol.coordinate.first, detail::to_simd_if<TypeParam>(1));
    SIMD_OR_SCALAR_EQ(aol.coordinate.second, detail::to_simd_if<TypeParam>(4));

    alignment_coordinate<TypeParam> coordinate_r{};
    coordinate_r.first = detail::to_simd_if<TypeParam>(2);
    coordinate_r.second = detail::to_simd_if<TypeParam>(3);
    detail::alignment_optimum aor{detail::to_simd_if<TypeParam>(10), coordinate_r};
    SIMD_OR_SCALAR_EQ(aor.score, detail::to_simd_if<TypeParam>(10));
    SIMD_OR_SCALAR_EQ(aor.coordinate.first, detail::to_simd_if<TypeParam>(2));
    SIMD_OR_SCALAR_EQ(aor.coordinate.second, detail::to_simd_if<TypeParam>(3));

    auto val = max(aol, aor);
    SIMD_OR_SCALAR_EQ(val.score, detail::to_simd_if<TypeParam>(10));
    SIMD_OR_SCALAR_EQ(val.coordinate.first, detail::to_simd_if<TypeParam>(2));
    SIMD_OR_SCALAR_EQ(val.coordinate.second, detail::to_simd_if<TypeParam>(3));

    if constexpr (Simd<TypeParam>)
    {
        if constexpr (simd_traits<TypeParam>::length >= 4)
            {
            alignment_coordinate<TypeParam> coordinate_l{};
            coordinate_l.first = iota<TypeParam>(0);
            coordinate_l.second = iota<TypeParam>(0);
            alignment_coordinate<TypeParam> coordinate_r{};
            coordinate_r.first = iota<TypeParam>(4);
            coordinate_r.second = iota<TypeParam>(4);

            detail::alignment_optimum<TypeParam> opt_l{};
            opt_l.score[0] = 0; opt_l.score[1] = 1; opt_l.score[2] = 0; opt_l.score[3] = 1;
            opt_l.coordinate = coordinate_l;

            detail::alignment_optimum<TypeParam> opt_r{};
            opt_r.score[0] = 1; opt_r.score[1] = 0; opt_r.score[2] = 1; opt_r.score[3] = 1;
            opt_r.coordinate = coordinate_r;

            TypeParam cmp_score{};
            TypeParam cmp_first{};
            TypeParam cmp_second{};

            cmp_score[0]  = 1; cmp_score[1]  = 1; cmp_score[2]  = 1; cmp_score[3]  = 1;
            cmp_first[0]  = 4; cmp_first[1]  = 1; cmp_first[2]  = 6; cmp_first[3]  = 3;
            cmp_second[0] = 4; cmp_second[1] = 1; cmp_second[2] = 6; cmp_second[3] = 3;

            auto res = max(opt_l, opt_r);

            SIMD_EQ(res.score, cmp_score);
            SIMD_EQ(res.coordinate.first, cmp_first);
            SIMD_EQ(res.coordinate.second, cmp_second);
        }
    }
}
