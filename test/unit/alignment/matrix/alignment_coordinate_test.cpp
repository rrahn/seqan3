// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <range/v3/view/iota.hpp>

#include <seqan3/alignment/detail/misc.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/core/simd/all.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

template <typename t>
struct advanceable_alignment_coordinate_test : public ::testing::Test
{};

template <typename t>
struct alignment_coordinate_test : public ::testing::Test
{};

using testing_types = ::testing::Types<int32_t, simd_type_t<int32_t>>;
TYPED_TEST_CASE(advanceable_alignment_coordinate_test, testing_types);
TYPED_TEST_CASE(alignment_coordinate_test, testing_types);

TEST(advanceable_alignment_coordinate, column_index_type)
{
    detail::column_index_type ci{1u};
    EXPECT_EQ(ci.get(), static_cast<size_t>(1));
}

TEST(advanceable_alignment_coordinate, row_index_type)
{
    detail::row_index_type ri{1u};
    EXPECT_EQ(ri.get(), static_cast<size_t>(1));
}

TYPED_TEST(advanceable_alignment_coordinate_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible<detail::advanceable_alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_constructible<detail::advanceable_alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_assignable<detail::advanceable_alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_constructible<detail::advanceable_alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_assignable<detail::advanceable_alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_destructible<detail::advanceable_alignment_coordinate<TypeParam>>::value);
}

TYPED_TEST(advanceable_alignment_coordinate_test, construction_with_different_state)
{
    detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>
        ro{detail::column_index_type{2u}, detail::row_index_type{3u}};

    detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none> no{ro};

    SIMD_OR_SCALAR_EQ(no.first,  detail::to_simd_if<TypeParam>(2));
    SIMD_OR_SCALAR_EQ(no.second, detail::to_simd_if<TypeParam>(3));
}

TYPED_TEST(advanceable_alignment_coordinate_test, type_deduction)
{
    detail::advanceable_alignment_coordinate def_co{};
    EXPECT_TRUE((std::is_same_v<decltype(def_co),
                    detail::advanceable_alignment_coordinate<size_t,
                                                             detail::advanceable_alignment_coordinate_state::none>>));

    detail::advanceable_alignment_coordinate co{detail::column_index_type{2u}, detail::row_index_type{3u}};
    EXPECT_TRUE((std::is_same_v<decltype(co),
                    detail::advanceable_alignment_coordinate<size_t,
                                                             detail::advanceable_alignment_coordinate_state::none>>));
}

TYPED_TEST(advanceable_alignment_coordinate_test, access)
{
    detail::advanceable_alignment_coordinate<TypeParam> def_co{};
    SIMD_OR_SCALAR_EQ(def_co.first,  detail::to_simd_if<TypeParam>(0));
    SIMD_OR_SCALAR_EQ(def_co.second, detail::to_simd_if<TypeParam>(0));

    detail::advanceable_alignment_coordinate<TypeParam> co{detail::column_index_type{2u}, detail::row_index_type{3u}};
    SIMD_OR_SCALAR_EQ(co.first,  detail::to_simd_if<TypeParam>(2));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(3));
}

TYPED_TEST(advanceable_alignment_coordinate_test, weakly_equality_comparable_concept)
{
    using not_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_TRUE(std::EqualityComparable<not_incrementable>);
    EXPECT_TRUE(std::EqualityComparable<row_incrementable>);
    EXPECT_TRUE(std::EqualityComparable<column_incrementable>);
}

TYPED_TEST(advanceable_alignment_coordinate_test, equality)
{
    using test_type =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{detail::column_index_type{10u}, detail::row_index_type{5u}};
    test_type t2{detail::column_index_type{5u}, detail::row_index_type{5u}};
    test_type t3{detail::column_index_type{10u}, detail::row_index_type{10u}};

    EXPECT_TRUE(t1 == t1);
    EXPECT_FALSE(t2 == t1);
    EXPECT_FALSE(t1 == t3);
    EXPECT_FALSE(t2 == t3);
}

TYPED_TEST(advanceable_alignment_coordinate_test, inequality)
{
    using test_type =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none>;

    test_type t1{detail::column_index_type{10u}, detail::row_index_type{5u}};
    test_type t2{detail::column_index_type{5u}, detail::row_index_type{5u}};
    test_type t3{detail::column_index_type{10u}, detail::row_index_type{10u}};

    EXPECT_FALSE(t1 != t1);
    EXPECT_TRUE(t2 != t1);
    EXPECT_TRUE(t1 != t3);
    EXPECT_TRUE(t2 != t3);
}

TYPED_TEST(advanceable_alignment_coordinate_test, incremental_concept)
{
    using not_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;
    using column_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    EXPECT_FALSE(std::WeaklyIncrementable<not_incrementable>);
    EXPECT_TRUE(std::WeaklyIncrementable<row_incrementable>);
    EXPECT_TRUE(std::WeaklyIncrementable<column_incrementable>);
}

TYPED_TEST(advanceable_alignment_coordinate_test, increment_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = ++co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(1u));
    auto co_tmp = co++;
    SIMD_OR_SCALAR_EQ(co_tmp.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co_tmp.second, detail::to_simd_if<TypeParam>(1u));
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(2u));
    co += 4;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(6u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, increment_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = ++co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(1u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
    auto co_tmp = co++;
    SIMD_OR_SCALAR_EQ(co_tmp.first, detail::to_simd_if<TypeParam>(1u));
    SIMD_OR_SCALAR_EQ(co_tmp.second, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(2u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
    co += 4;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(6u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, decrement_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co += 4;
    auto co_tmp = co--;
    SIMD_OR_SCALAR_EQ(co_tmp.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co_tmp.second, detail::to_simd_if<TypeParam>(4u));
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(3u));

    co = --co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(2u));

    co -= 2;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, decrement_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co += 4;
    auto co_tmp = co--;
    SIMD_OR_SCALAR_EQ(co_tmp.first, detail::to_simd_if<TypeParam>(4u));
    SIMD_OR_SCALAR_EQ(co_tmp.second, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(3u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));

    co = --co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(2u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));

    co -= 2;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, advance_row)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};

    co = co + 4;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(4u));

    co = 4 + co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(0u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(8u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, advance_col)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co{detail::column_index_type{0u}, detail::row_index_type{0u}};
    co = co + 4;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(4u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));

    co = 4 + co;
    SIMD_OR_SCALAR_EQ(co.first, detail::to_simd_if<TypeParam>(8u));
    SIMD_OR_SCALAR_EQ(co.second, detail::to_simd_if<TypeParam>(0u));
}

TYPED_TEST(advanceable_alignment_coordinate_test, iota_column_index)
{
    using col_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    col_incrementable co_begin{detail::column_index_type{0u}, detail::row_index_type{0u}};
    col_incrementable co_end{detail::column_index_type{5u}, detail::row_index_type{0u}};
    auto v = std::view::iota(co_begin, co_end);

    EXPECT_TRUE((std::Same<decltype(v.begin()), decltype(v.end())>));
    SIMD_OR_SCALAR_EQ((*(--v.end())).first, detail::to_simd_if<TypeParam>(4u));

    size_t test = 0u;
    for (auto coordinate : v)
        SIMD_OR_SCALAR_EQ(coordinate.first, detail::to_simd_if<TypeParam>(test++));
}

TYPED_TEST(advanceable_alignment_coordinate_test, iota_row_index)
{
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;

    row_incrementable co_begin{detail::column_index_type{0u}, detail::row_index_type{0u}};
    row_incrementable co_end{detail::column_index_type{0u}, detail::row_index_type{5u}};
    auto v = std::view::iota(co_begin, co_end);

    EXPECT_TRUE((std::Same<decltype(v.begin()), decltype(v.end())>));
    SIMD_OR_SCALAR_EQ((*(--v.end())).second, detail::to_simd_if<TypeParam>(4u));

    size_t test = 0u;
    for (auto coordinate : v)
        SIMD_OR_SCALAR_EQ(coordinate.second, detail::to_simd_if<TypeParam>(test++));
}

// TODO: Need to be reenabled once the issue with the debug_stream is fixed.
// TEST(advanceable_alignment_coordinate, debug_stream)
// {
//     using not_incrementable =
//         detail::advanceable_alignment_coordinate<size_t, detail::advanceable_alignment_coordinate_state::none>;
//     using row_incrementable =
//         detail::advanceable_alignment_coordinate<size_t, detail::advanceable_alignment_coordinate_state::row>;
//     using col_incrementable =
//         detail::advanceable_alignment_coordinate<size_t, detail::advanceable_alignment_coordinate_state::column>;
//
//     not_incrementable co_not{detail::column_index_type{10u}, detail::row_index_type{5u}};
//     col_incrementable co_col{detail::column_index_type{10u}, detail::row_index_type{5u}};
//     row_incrementable co_row{detail::column_index_type{10u}, detail::row_index_type{5u}};
//
//     std::stringstream sstream{};
//     debug_stream_type dstream{sstream};
//     dstream << co_not;
//     dstream << co_col;
//     dstream << co_row;
//     EXPECT_EQ(sstream.str(), "(10,5)(10,5)(10,5)");
//
//     EXPECT_EQ(co_not, co_not);
//     EXPECT_EQ(co_col, co_col);
//     EXPECT_EQ(co_row, co_row);
// }

TYPED_TEST(alignment_coordinate_test, basic)
{
    EXPECT_TRUE(std::is_default_constructible<alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_constructible<alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_assignable<alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_constructible<alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_assignable<alignment_coordinate<TypeParam>>::value);
    EXPECT_TRUE(std::is_destructible<alignment_coordinate<TypeParam>>::value);

    using not_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::row>;
    using col_incrementable =
        detail::advanceable_alignment_coordinate<TypeParam, detail::advanceable_alignment_coordinate_state::column>;

    not_incrementable co_not{detail::column_index_type{10u}, detail::row_index_type{5u}};
    col_incrementable co_col{detail::column_index_type{10u}, detail::row_index_type{5u}};
    row_incrementable co_row{detail::column_index_type{10u}, detail::row_index_type{5u}};

    alignment_coordinate<TypeParam> test1{co_not};
    SIMD_OR_SCALAR_EQ(test1.first, detail::to_simd_if<TypeParam>(10u));
    SIMD_OR_SCALAR_EQ(test1.second, detail::to_simd_if<TypeParam>(5u));

    alignment_coordinate<TypeParam> test2{co_col};
    SIMD_OR_SCALAR_EQ(test2.first, detail::to_simd_if<TypeParam>(10u));
    SIMD_OR_SCALAR_EQ(test2.second, detail::to_simd_if<TypeParam>(5u));

    alignment_coordinate<TypeParam> test3{co_row};
    SIMD_OR_SCALAR_EQ(test3.first, detail::to_simd_if<TypeParam>(10u));
    SIMD_OR_SCALAR_EQ(test3.second, detail::to_simd_if<TypeParam>(5u));

    alignment_coordinate<TypeParam> test4{detail::column_index_type{10u}, detail::row_index_type{5u}};
    SIMD_OR_SCALAR_EQ(test4.first, detail::to_simd_if<TypeParam>(10u));
    SIMD_OR_SCALAR_EQ(test4.second, detail::to_simd_if<TypeParam>(5u));
}

// TODO Enable once debug stream problem is fixed.
// TEST(alignment_coordinate, debug_stream)
// {
//     alignment_coordinate co_align{detail::column_index_type{10u}, detail::row_index_type{5u}};
//
//     std::stringstream sstream{};
//     debug_stream_type dstream{sstream};
//     dstream << co_align;
//     EXPECT_EQ(sstream.str(), "(10,5)");
//
//     EXPECT_EQ(co_align, co_align);
// }
