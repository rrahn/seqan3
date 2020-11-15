// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/alignment/matrix/detail/trace_matrix_simd_adapter_iterator.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/test/simd_utility.hpp>

#include "../../../range/iterator_test_template.hpp"

using trace_t = seqan3::simd::simd_type_t<uint8_t>; // trace_directions
using allocator_t = seqan3::aligned_allocator<trace_t, sizeof(trace_t)>;
using matrix_t = seqan3::detail::two_dimensional_matrix<trace_t,
                                                        allocator_t,
                                                        seqan3::detail::matrix_major_order::column>;

// ----------------------------------------------------------------------------
// two dimensional simd matrix iterator
// ----------------------------------------------------------------------------

using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<trace_t, allocator_t> expected_range{20, trace_t{}};
    matrix_t test_range;

    void SetUp()
    {
        test_range.resize(seqan3::detail::number_rows{5}, seqan3::detail::number_cols{4});
    }

    static void expect_eq(trace_t const actual_cell, trace_t const expected_cell)
    {
        SIMD_EQ(actual_cell, expected_cell);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(two_dimensional_simd_matrix_iterator_test, iterator_fixture, matrix_iterator_t, );

// ----------------------------------------------------------------------------
// two dimensional simd matrix iterator adapter
// ----------------------------------------------------------------------------

using adapter_iterator_t = seqan3::detail::trace_matrix_simd_adapter_iterator<matrix_iterator_t>;

template <>
struct iterator_fixture<adapter_iterator_t> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    using test_range_t = std::ranges::subrange<adapter_iterator_t, adapter_iterator_t>;

    std::vector<seqan3::detail::trace_directions> expected_range{20, seqan3::detail::trace_directions::none};
    matrix_t test_range_data;
    test_range_t test_range{};

    void SetUp()
    {
        test_range_data.resize(seqan3::detail::number_rows{5}, seqan3::detail::number_cols{4});
        test_range = test_range_t{adapter_iterator_t{test_range_data.begin(), 0},
                                  adapter_iterator_t{test_range_data.end(), 0}};
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(two_dimensional_simd_matrix_iterator_adapter_test,
                               iterator_fixture,
                               adapter_iterator_t, );

TEST(two_dimensional_simd_matrix_iterator_adapter_test, concept)
{
    EXPECT_TRUE(seqan3::detail::two_dimensional_matrix_iterator<adapter_iterator_t>);
}

TEST(two_dimensional_simd_matrix_iterator_adapter_test, trace_path)
{
    using dir_t = seqan3::detail::trace_directions;

    matrix_t matrix{};
    matrix.resize(seqan3::detail::number_rows{3}, seqan3::detail::number_cols{4});
    auto trace_matrix_it = matrix.begin();

    // Initialise column 0
    *trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::none));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::up_open));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::up));

    // Initialise column 1
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::left_open));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::diagonal));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::up_open));

    // Initialise column 2
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::left));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::diagonal));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::left_open));

    // Initialise column 3
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::left));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::up_open));
    *++trace_matrix_it = seqan3::simd::fill<trace_t>(static_cast<uint8_t>(dir_t::left));

    EXPECT_TRUE(++trace_matrix_it == matrix.end());

    trace_matrix_it = matrix.begin() + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                     seqan3::detail::column_index_type{3}};

    using trace_iterator_t = seqan3::detail::trace_iterator<adapter_iterator_t>;
    using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

    adapter_iterator_t adapter_iter{trace_matrix_it, 0};
    path_t trace_path{trace_iterator_t{adapter_iter}, std::default_sentinel};

    auto trace_path_it = trace_path.begin();
    EXPECT_EQ(*trace_path_it, dir_t::left);
    EXPECT_EQ(*++trace_path_it, dir_t::left);
    EXPECT_EQ(*++trace_path_it, dir_t::up);
    EXPECT_EQ(*++trace_path_it, dir_t::diagonal);
    EXPECT_EQ(*++trace_path_it, dir_t::none);
    EXPECT_TRUE(trace_path_it == trace_path.end());
}
