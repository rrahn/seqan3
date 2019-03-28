// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/gap_scheme_simd_wrapper.hpp>
#include <seqan3/test/simd_utility.hpp>

#if SEQAN3_WITH_CEREAL
#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>

#include <seqan3/test/tmp_filename.hpp>
#endif // SEQAN3_WITH_CEREAL

using namespace seqan3;
using namespace seqan3::detail;

TEST(gap_scheme_simd_wrapper, construction)
{
    using simd_t = typename simd_type<int32_t>::type;
    using scheme_t = gap_scheme_simd_wrapper<simd_t>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_destructible_v<scheme_t>);
    EXPECT_TRUE((std::is_constructible_v<scheme_t, gap_scheme<int32_t>>));
}

TEST(gap_scheme_simd_wrapper, deduction)
{
    {
        using simd_t = typename simd_type<int8_t>::type;
        gap_scheme_simd_wrapper scheme{gap_scheme{gap_score{-1}, gap_open_score{-10}}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme_simd_wrapper<simd_t>>));
    }

    { // TODO: Enable simd for floating point numbers
        // using simd_t = typename simd_type<float>::type;
        // gap_scheme_simd_wrapper scheme{gap_scheme{gap_score{-1.0}, gap_open_score{-10.0}}};
        // EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme_simd_wrapper<simd_t>>));
    }
}

TEST(gap_scheme_simd_wrapper, set_gap_scheme)
{
    {  // exception: value gap to causes overflow
        using simd_t = typename simd_type<int8_t>::type;
        using scheme_t = gap_scheme_simd_wrapper<simd_t>;
        EXPECT_THROW((scheme_t{gap_scheme<int16_t>{gap_score{-300}, gap_open_score{-400}}}), std::invalid_argument);
    }

    { // exception: value gap open causes overflow
        using simd_t = typename simd_type<int8_t>::type;
        using scheme_t = gap_scheme_simd_wrapper<simd_t>;
        EXPECT_THROW((scheme_t{gap_scheme<int16_t>{gap_score{-255}, gap_open_score{-300}}}), std::invalid_argument);
    }

    { // set the gap_scheme
        using simd_t = typename simd_type<int32_t>::type;
        gap_scheme_simd_wrapper<simd_t> scheme{gap_scheme{gap_score{-1}, gap_open_score{-10}}};

        SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-1));
        SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(-10));
    }
}

TEST(simd_gap_scheme_wrapper, member_types)
{
    using simd_t = simd_type_t<int32_t>;
    gap_scheme_simd_wrapper<simd_t> scheme{};

    using score_t = typename gap_scheme_simd_wrapper<simd_t>::score_type;
    EXPECT_TRUE((std::is_same_v<score_t, simd_t>));
}

TEST(simd_gap_scheme_wrapper, get_gap_score)
{
    using simd_t = simd_type_t<int32_t>;

    gap_scheme_simd_wrapper<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-1));
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_score())>));
}


TEST(simd_gap_scheme_wrapper, get_gap_open_score)
{
    using simd_t = simd_type_t<int32_t>;

    gap_scheme_simd_wrapper<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(0));
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_open_score())>));
}

// TODO Currently not working.
// #if SEQAN3_WITH_CEREAL
// template <typename in_archive_t, typename out_archive_t, typename TypeParam>
// void do_serialisation(TypeParam const l)
// {
//     // This makes sure the file is also deleted if an exception is thrown in one of the tests below
//     // Generate unique file name.
//     test::tmp_filename filename{"simd_gap_scheme_cereal_test"};
//
//     {
//         std::ofstream os{filename.get_path(), std::ios::binary};
//         out_archive_t oarchive{os};
//         oarchive(l);
//     }
//
//     {
//         TypeParam in_l{};
//         std::ifstream is{filename.get_path(), std::ios::binary};
//         in_archive_t iarchive{is};
//         iarchive(in_l);
//         EXPECT_EQ(l, in_l);
//     }
// }
//
// TEST(simd_gap_scheme_wrapper, serialisation)
// {
//     using simd_t = simd_type_t<int32_t>;
//
//     gap_scheme_simd_wrapper<simd_t> scheme1;
//     scheme1.set_linear(gap_score{-3});
//
//     do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
//     do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
//     do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
//     do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
//
//     scheme1.set_affine(gap_score{-3}, gap_open_score{-6});
//
//     do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
//     do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
//     do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
//     do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
// }
// #endif // SEQAN3_WITH_CEREAL
