// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/std/span>
#include <seqan3/test/simd_utility.hpp>

template <typename t>
class simd_algorithm_t : public ::testing::Test
{
public:
};

using test_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t>;

TYPED_TEST_CASE(simd_algorithm_t, test_types);

using namespace seqan3;

TEST(simd_algorithm, fill)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = 4;

    constexpr simd_type result = fill<simd_type>(4);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = i;

    constexpr simd_type result = iota<simd_type>(0);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, transform_to_soa)
{
    using simd_t = simd_type_t<int32_t, 4>;

    EXPECT_EQ(simd_traits<simd_t>::length, 4);

    auto seq1 = "ATGCAAAAATA"_dna4;
    auto seq2 = "CATGCCCCCGC"_dna4;
    auto seq3 = "GCATGGGGGGC"_dna4;
    auto seq4 = "TGCATTTTTTA"_dna4;

    std::vector data{std::span{std::ranges::data(seq1), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq2), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq3), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq4), simd_traits<simd_t>::length}};

    std::vector<simd_t, aligned_allocator<simd_t, simd_traits<simd_t>::max_length>> out_vec;

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[0], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[1], simd_t{3, 0, 1, 2});
    SIMD_EQ(out_vec[2], simd_t{2, 3, 0, 1});
    SIMD_EQ(out_vec[3], simd_t{1, 2, 3, 0});

    data[0] = std::span{std::ranges::data(seq1) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[1] = std::span{std::ranges::data(seq2) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[2] = std::span{std::ranges::data(seq3) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[3] = std::span{std::ranges::data(seq4) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[4], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[5], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[6], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[7], simd_t{0, 1, 2, 3});

    data[0] = std::span{std::ranges::data(seq1) + 2 * simd_traits<simd_t>::length, 2};
    data[1] = std::span{std::ranges::data(seq2) + 2 * simd_traits<simd_t>::length, 2};
    data[2] = std::span{std::ranges::data(seq3) + 2 * simd_traits<simd_t>::length, 2};
    data[3] = std::span{std::ranges::data(seq4) + 2 * simd_traits<simd_t>::length, 2};

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[8],  simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[9],  simd_t{3, 2, 2, 3});
    SIMD_EQ(out_vec[10], simd_t{0, 1, 1, 0});
}

TYPED_TEST(simd_algorithm_t, unpack_hi)
{
    using simd_t = typename simd_type<TypeParam>::type;

    // Nothing to test if the simd type is just scalar.
    // TODO: report gcc9.1 fails for intrinsic _mm512_permutex2var_epi8 on mac: g++-9 (Homebrew GCC 9.1.0) 9.1.0
    if constexpr (simd_traits<simd_t>::length == 1 || simd_traits<simd_t>::length == 64)
    {
        return;
    }
    else
    {
        simd_t lhs = iota<simd_t>(1);
        simd_t rhs = iota<simd_t>(1);

        simd_t res = unpack_hi(lhs, rhs);

        simd_t cmp{};

        uint32_t j = simd_traits<simd_t>::length / 2;
        for (uint32_t i = 0; i < simd_traits<simd_t>::length; i += 2, ++j)
        {
            cmp[i]     = lhs[j];
            cmp[i + 1] = rhs[j];
        }

        SIMD_EQ(res, cmp);
    }
}

TYPED_TEST(simd_algorithm_t, unpack_lo)
{
    using simd_t = typename simd_type<TypeParam>::type;

    // Nothing to test if the simd type is just scalar.
    // TODO: report gcc9.1 fails for intrinsic _mm512_permutex2var_epi8 on mac: g++-9 (Homebrew GCC 9.1.0) 9.1.0
    if constexpr (simd_traits<simd_t>::length == 1 || simd_traits<simd_t>::length == 64)
    {
        return;
    }
    else
    {

        simd_t lhs = iota<simd_t>(1);
        simd_t rhs = iota<simd_t>(1);

        simd_t res = unpack_lo(lhs, rhs);

        simd_t cmp{};

        for (uint32_t i = 0, j = 0; i < simd_traits<simd_t>::length; i += 2, ++j)
        {
            cmp[i]     = lhs[j];
            cmp[i + 1] = rhs[j];
        }

        SIMD_EQ(res, cmp);
    }
}

TEST(simd_algorithm, transpose)
{
    using simd_t = simd_type_t<uint8_t>;

    if constexpr (simd_traits<simd_t>::length > 1)
    {
        std::array<simd_t, simd_traits<simd_t>::length> matrix;

        for (size_t i = 0; i < matrix.size(); ++i)
            matrix[i] = simd::iota<simd_t>(0);

        simd::transpose(matrix);

        for (size_t i = 0; i < matrix.size(); ++i)
            SIMD_EQ(matrix[i], simd::fill<simd_t>(i));
    }
}

template <typename t>
class simd_algorithm_upcast : public ::testing::Test
{
public:

    static constexpr auto target_list_signed()
    {
        if constexpr (sizeof(t) == 1)
            return type_list_signed_8{};
        else if constexpr (sizeof(t) == 2)
            return type_list_signed_16{};
        else
            return type_list_signed_32{};
    }

    static constexpr auto target_list_unsigned()
    {
        if constexpr (sizeof(t) == 1)
            return type_list_unsigned_8{};
        else if constexpr (sizeof(t) == 2)
            return type_list_unsigned_16{};
        else
            return type_list_unsigned_32{};
    }

    using type_list_signed_8  = type_list<int16_t, int32_t, int64_t>;
    using type_list_signed_16 = type_list<int32_t, int64_t>;
    using type_list_signed_32 = type_list<int64_t>;

    using type_list_unsigned_8  = type_list<uint16_t, uint32_t, uint64_t>;
    using type_list_unsigned_16 = type_list<uint32_t, uint64_t>;
    using type_list_unsigned_32 = type_list<uint64_t>;
};

using upcast_test_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t>;
TYPED_TEST_CASE(simd_algorithm_upcast, upcast_test_types);

TYPED_TEST(simd_algorithm_upcast, signed)
{
    using list = decltype(TestFixture::target_list_signed());

    detail::for_each_type([this] (auto type_id)
    {
        using target_type = typename decltype(type_id)::type;
        using src_simd_t = simd_type_t<TypeParam>;
        using tar_simd_t = simd_type_t<target_type>;

        src_simd_t s = fill<src_simd_t>(-10);
        tar_simd_t t = upcast_auto<tar_simd_t>(s);

        for (size_t i  = 0; i < simd_traits<tar_simd_t>::length; ++i)
            EXPECT_EQ(t[i], static_cast<target_type>(static_cast<TypeParam>(-10)));
    }, list{});
}

TYPED_TEST(simd_algorithm_upcast, unsigned)
{
    using list = decltype(TestFixture::target_list_unsigned());

    detail::for_each_type([this] (auto type_id)
    {
        using target_type = typename decltype(type_id)::type;
        using src_simd_t = simd_type_t<TypeParam>;
        using tar_simd_t = simd_type_t<target_type>;

        src_simd_t s = fill<src_simd_t>(-10);
        tar_simd_t t = upcast_auto<tar_simd_t>(s);

        for (size_t i  = 0; i < simd_traits<tar_simd_t>::length; ++i)
            EXPECT_EQ(t[i], static_cast<target_type>(static_cast<TypeParam>(-10)));
    }, list{});
}
