// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/view_iota_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/test/simd_utility.hpp>

#include "../../range/iterator_test_template.hpp"

template <typename simd_t>
struct iota_simd_test : public ::testing::Test
{};

template <typename simd_t>
struct iterator_fixture<iota_simd_test<simd_t>> : public iota_simd_test<simd_t>
{
    using iota_simd_view_type = seqan3::views::iota_simd<simd_t>;;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<simd_t, seqan3::aligned_allocator<simd_t, alignof(simd_t)>> expected_range{};
    iota_simd_view_type test_range{};

    void SetUp()
    {
        using scalar_t = typename seqan3::simd::simd_traits<simd_t>::scalar_type;
        for (scalar_t i = 0; i < static_cast<scalar_t>(255); ++i)
            expected_range.push_back(seqan3::simd::fill<simd_t>(i));

        test_range = iota_simd_view_type{static_cast<scalar_t>(0), static_cast<scalar_t>(255)};
    }

    template <typename actual_simd_t, typename expected_simd_t>
    static void expect_eq(actual_simd_t && actual_index, expected_simd_t && expected_index)
    {
        SIMD_EQ(actual_index, expected_index);
    }
};

using testing_types = ::testing::Types<iota_simd_test<seqan3::simd::simd_type_t<uint8_t>>,
                                       iota_simd_test<seqan3::simd::simd_type_t<uint16_t>>,
                                       iota_simd_test<seqan3::simd::simd_type_t<int32_t>>,
                                       iota_simd_test<seqan3::simd::simd_type_t<int64_t>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(view_iota_simd_iterator, iterator_fixture, testing_types, );

TEST(iota_simd_test, concept_checks)
{
    using simd_t = seqan3::simd::simd_type_t<uint32_t>;
    using iota_simd_view_t = seqan3::views::iota_simd<simd_t>;

    EXPECT_TRUE(std::ranges::forward_range<iota_simd_view_t>);
    EXPECT_TRUE(std::ranges::common_range<iota_simd_view_t>);
    // EXPECT_TRUE(std::ranges::borrowed_range<iota_simd_view_t>);
    EXPECT_TRUE(std::ranges::sized_range<iota_simd_view_t>);
}

TEST(iota_simd_test, size)
{
    using simd_t = seqan3::simd::simd_type_t<uint32_t>;
    seqan3::views::iota_simd<simd_t> test_view{0, 10};

    EXPECT_EQ(test_view.size(), 10u);
}
