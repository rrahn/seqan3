// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_identity.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_storage_two_dimensional.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix.hpp>
#include <seqan3/test/expect_range_eq.hpp>

struct alignment_matrix_test : public ::testing::Test
{
    using element_t = seqan3::detail::alignment_matrix_element_identity<int32_t>;
    using storage_t = seqan3::detail::alignment_matrix_storage_two_dimensional<element_t>;
    using matrix_t = seqan3::detail::alignment_matrix<storage_t>;

    seqan3::detail::row_index_type<size_t> rdim{4};
    seqan3::detail::column_index_type<size_t> cdim{5};
    storage_t s{};
    matrix_t m{};
};

TEST_F(alignment_matrix_test, range_concept)
{
    EXPECT_TRUE(std::ranges::input_range<matrix_t>);
}

TEST_F(alignment_matrix_test, storage_get)
{
    EXPECT_EQ(m.rows(), 0u);
    EXPECT_EQ(m.cols(), 0u);

    m.storage().resize(rdim, cdim);

    EXPECT_EQ(m.rows(), rdim.get());
    EXPECT_EQ(m.cols(), cdim.get());
}

TEST_F(alignment_matrix_test, storage_set)
{
    storage_t s;
    s.resize(rdim, cdim);
    m.storage(std::as_const(s));

    EXPECT_EQ(m.rows(), rdim.get());
    EXPECT_EQ(m.cols(), cdim.get());

    m.storage().clear();

    EXPECT_EQ(m.rows(), 0u);
    EXPECT_EQ(m.cols(), 0u);

    m.storage(std::move(s));

    EXPECT_EQ(m.rows(), rdim.get());
    EXPECT_EQ(m.cols(), cdim.get());
}

TEST_F(alignment_matrix_test, iterate)
{
    m.storage().resize(rdim, cdim, 10);

    for (auto column : m)
        EXPECT_RANGE_EQ(column, std::vector(rdim.get(), 10));
}
