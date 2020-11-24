// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_identity.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_storage_concept.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_storage_sparse.hpp>
#include <seqan3/test/expect_range_eq.hpp>

struct alignment_matrix_storage_sparse_test : public ::testing::Test
{
    using element_t = seqan3::detail::alignment_matrix_element_identity<int32_t>;
    using storage_t = seqan3::detail::alignment_matrix_storage_sparse<element_t>;

    seqan3::detail::row_index_type<size_t> rdim{4};
    seqan3::detail::column_index_type<size_t> cdim{5};
    storage_t s{};
};

TEST_F(alignment_matrix_storage_sparse_test, matrix_storage_concept)
{
    EXPECT_TRUE(seqan3::detail::alignment_matrix_storage<storage_t>);
}

TEST_F(alignment_matrix_storage_sparse_test, cols)
{
    s.resize(rdim, cdim);
    EXPECT_EQ(s.cols(), cdim.get());
}

TEST_F(alignment_matrix_storage_sparse_test, rows)
{
    s.resize(rdim, cdim);
    EXPECT_EQ(s.rows(), rdim.get());
}

TEST_F(alignment_matrix_storage_sparse_test, column_at)
{
    s.resize(rdim, cdim, 10);

    EXPECT_RANGE_EQ(s.column_at(0), std::vector(rdim.get(), 10));
    EXPECT_RANGE_EQ(s.column_at(1), std::vector(rdim.get(), 10));
    EXPECT_RANGE_EQ(s.column_at(2), std::vector(rdim.get(), 10));
    EXPECT_RANGE_EQ(s.column_at(3), std::vector(rdim.get(), 10));
    EXPECT_RANGE_EQ(s.column_at(4), std::vector(rdim.get(), 10));
}

TEST_F(alignment_matrix_storage_sparse_test, clear)
{
    s.resize(rdim, cdim);
    s.clear();
    EXPECT_EQ(s.cols(), 0U);
    EXPECT_EQ(s.rows(), 0U);
}
