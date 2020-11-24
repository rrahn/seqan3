// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_concept.hpp>
#include <seqan3/alignment/matrix/detail/alignment_matrix_element_identity.hpp>

struct alignment_matrix_element_identity_test : public ::testing::Test
{
    using element_t = seqan3::detail::alignment_matrix_element_identity<int32_t>;

    element_t e{};
};

TEST_F(alignment_matrix_element_identity_test, matrix_element_concept)
{
    EXPECT_TRUE(seqan3::detail::alignment_matrix_element<element_t>);
}

TEST_F(alignment_matrix_element_identity_test, value_type)
{
    // Test that the value type is the passed template parameter.
    EXPECT_TRUE((std::same_as<typename element_t::value_type, int32_t>));
}

TEST_F(alignment_matrix_element_identity_test, storage_value_type)
{
    // In the identity matrix element the storage value type is identical to the value type.
    using value_t =  typename element_t::value_type;
    EXPECT_TRUE((std::same_as<typename element_t::storage_value_type, value_t>));
}

TEST_F(alignment_matrix_element_identity_test, element_type)
{
    // In the identity matrix element the element type is identical to the value type.
    using value_t =  typename element_t::value_type;
    EXPECT_TRUE((std::same_as<typename element_t::element_type, value_t>));
}

TEST_F(alignment_matrix_element_identity_test, element_reference)
{
    // In the identity matrix element the element type is identical to the lvalue reference of the value type.
    using value_t =  typename element_t::value_type;
    EXPECT_TRUE((std::same_as<typename element_t::element_reference, std::add_lvalue_reference_t<value_t>>));
}

TEST_F(alignment_matrix_element_identity_test, make_element)
{
    // Make alignment is called by the matrix storage iterator to return the element_reference as defined
    // by the given alignment matrix element.
    using storage_value_t = typename element_t::storage_value_type;
    storage_value_t storage_value{1};
    typename element_t::element_reference ref = e.make_element(storage_value);
    EXPECT_EQ(ref, 1);
    ref = -10;
    EXPECT_EQ(ref, -10);
}

TEST_F(alignment_matrix_element_identity_test, initialise)
{
    // The alignment matrix storage might need to be initialised with a particular value.
    // To work properly, the initialisation is delegated to the matrix element type.
    using storage_value_t = typename element_t::storage_value_type;
    storage_value_t init_value = e.initialise(10);
    EXPECT_EQ(init_value, 10);
}
