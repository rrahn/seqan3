// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_concept.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_element_affine_cell.hpp>

struct alignment_matrix_element_affine_cell_test : public ::testing::Test
{
    using element_t = seqan3::detail::alignment_matrix_element_affine_cell<int32_t>;

    element_t e{};
};

TEST_F(alignment_matrix_element_affine_cell_test, matrix_element_concept)
{
    EXPECT_TRUE(seqan3::detail::alignment_matrix_element<element_t>);
}

TEST_F(alignment_matrix_element_affine_cell_test, value_type)
{
    // Test that the value type is the passed template parameter.
    EXPECT_TRUE((std::same_as<typename element_t::value_type, int32_t>));
}

TEST_F(alignment_matrix_element_affine_cell_test, storage_value_type)
{
    // Test that the value type for the matrix storage is usable for the affine alignment recursion.
    using value_t =  typename element_t::value_type;
    EXPECT_TRUE((std::same_as<typename element_t::storage_value_type, std::pair<value_t, value_t>>));
}

TEST_F(alignment_matrix_element_affine_cell_test, element_type)
{
    // The element type is the actual type of the alignment matrix.
    using value_t =  typename element_t::value_type;
    using affine_cell_t = seqan3::detail::affine_cell<std::tuple<value_t, value_t, value_t>>;
    EXPECT_TRUE((std::same_as<typename element_t::element_type, affine_cell_t>));
}

TEST_F(alignment_matrix_element_affine_cell_test, element_reference)
{
    // The element reference is the proxy returned by the iterator over the matrix elements.
    using value_t =  typename element_t::value_type;
    using affine_cell_proxy_t = seqan3::detail::affine_cell<std::tuple<value_t &, value_t &, value_t &>>;
    EXPECT_TRUE((std::same_as<typename element_t::element_reference, affine_cell_proxy_t>));
}

TEST_F(alignment_matrix_element_affine_cell_test, make_element)
{
    // Make alignment is called by the matrix storage iterator to return the element_reference as defined
    // by the given alignment matrix element
    using storage_value_t = typename element_t::storage_value_type;
    storage_value_t storage_value{1, 5};
    typename element_t::element_reference affine_cell = e.make_element(storage_value);
    EXPECT_EQ(affine_cell.best_score(), 1);
    EXPECT_EQ(affine_cell.horizontal_score(), 5);
    EXPECT_EQ(affine_cell.vertical_score(), 0);

    affine_cell.best_score() = -10;
    affine_cell.horizontal_score() = 0;
    affine_cell.vertical_score() = 100;

    EXPECT_EQ(affine_cell.best_score(), -10);
    EXPECT_EQ(affine_cell.horizontal_score(), 0);
    EXPECT_EQ(affine_cell.vertical_score(), 100);

    EXPECT_EQ(storage_value.first, -10);
    EXPECT_EQ(storage_value.second, 0);
}

TEST_F(alignment_matrix_element_affine_cell_test, initialise)
{
    // The alignment matrix storage might need to be initialised with a particular value.
    // To work properly, the initialisation is delegated to the matrix element type.
    using storage_value_t = typename element_t::storage_value_type;
    storage_value_t init_value = e.initialise(10);
    EXPECT_EQ(init_value.first, 10);
    EXPECT_EQ(init_value.second, 10);
    EXPECT_EQ(e.make_element(init_value).vertical_score(), 10);
}
