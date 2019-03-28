// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/simd/simd_debug_stream.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/view/to_simd_chunk.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

template <typename t>
class view_to_simd_chunk_t : public ::testing::Test
{
public:
    using container_t = std::tuple_element_t<0, t>;
    using simd_t      = std::tuple_element_t<1, t>;

    void SetUp()
    {
        data.resize(simd_traits<simd_t>::length);
        for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        {
            // Generate sequences that end on different boundaries
            size_t l = simd_traits<simd_t>::length * 64 - (i * simd_traits<simd_t>::length) - i;
            std::ranges::copy(test::generate_sequence<value_type_t<container_t>>(l), std::back_inserter(data[i]));
        }

        cmp_data.resize(simd_traits<simd_t>::length * 64);  // longest sequence in set.

        for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        {
            for (size_t j = 0; j < data[i].size(); ++j)
                cmp_data[j][i] = seqan3::to_rank(data[i][j]);
        }
    }

    std::vector<container_t> data;
    std::vector<simd_t, aligned_allocator<simd_t, alignof(simd_t)>> cmp_data;
};

using test_types = ::testing::Types<std::tuple<std::vector<dna4>, simd_type_t<int8_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int16_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int32_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int64_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint8_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint16_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint32_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint64_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int8_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int16_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int32_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int64_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint8_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint16_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint32_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint64_t>>
                                   >;
TYPED_TEST_CASE(view_to_simd_chunk_t, test_types);

TEST(view_to_simd_chunk, concept)
{
    using cmp_type = std::vector<dna4_vector>;
    using test_type = detail::view_to_simd_chunk<std::ranges::all_view<cmp_type &>, simd_type_t<int8_t>>;

    using iter_t = decltype(std::ranges::begin(std::declval<test_type &>()));
    EXPECT_TRUE(std::InputIterator<iter_t>);

    EXPECT_EQ(std::ranges::InputRange<cmp_type>, std::ranges::InputRange<test_type>);
    EXPECT_NE(std::ranges::ForwardRange<cmp_type>, std::ranges::ForwardRange<test_type>);
    EXPECT_NE(std::ranges::BidirectionalRange<cmp_type>, std::ranges::BidirectionalRange<test_type>);
    EXPECT_NE(std::ranges::RandomAccessRange<cmp_type>, std::ranges::RandomAccessRange<test_type>);
    EXPECT_NE(std::ranges::RandomAccessRange<cmp_type>, std::ranges::RandomAccessRange<test_type>);

    EXPECT_EQ(std::ranges::Range<cmp_type>, std::ranges::Range<test_type>);
    EXPECT_NE(std::ranges::View<cmp_type>, std::ranges::View<test_type>);
    EXPECT_EQ(std::ranges::SizedRange<cmp_type>, std::ranges::SizedRange<test_type>);
    EXPECT_NE(std::ranges::CommonRange<cmp_type>, std::ranges::CommonRange<test_type>);
    EXPECT_NE(ConstIterableRange<cmp_type>, ConstIterableRange<test_type>);
    EXPECT_NE((std::ranges::OutputRange<cmp_type, dna4_vector>), (std::ranges::OutputRange<test_type, dna4_vector>));
}

TEST(view_to_simd_chunk, iter_concept)
{
    using cmp_type = std::vector<dna4_vector>;
    using test_type = detail::view_to_simd_chunk<std::ranges::all_view<cmp_type &>, simd_type_t<int8_t>>;
    using iter_t = std::ranges::iterator_t<test_type>;
    using sent_t = std::ranges::sentinel_t<test_type>;

    EXPECT_TRUE(std::Iterator<iter_t>);
    EXPECT_TRUE(std::InputIterator<iter_t>);
    EXPECT_FALSE(std::ForwardIterator<iter_t>);
    EXPECT_FALSE(std::BidirectionalIterator<iter_t>);
    EXPECT_FALSE(std::RandomAccessIterator<iter_t>);
    EXPECT_FALSE((std::OutputIterator<iter_t, decltype(*std::declval<iter_t &>())>));
    EXPECT_TRUE((std::Sentinel<sent_t, iter_t>));
}

TYPED_TEST(view_to_simd_chunk_t, begin)
{
    auto v = this->data | view::to_simd_chunk<typename TestFixture::simd_t>;

    auto it = v.begin();
    auto chunk = *it;

    for (size_t i = 0; i < chunk.size(); ++i)
        SIMD_EQ(chunk[i], this->cmp_data[i]);
}

TYPED_TEST(view_to_simd_chunk_t, end)
{

}

TYPED_TEST(view_to_simd_chunk_t, size)
{

}

TYPED_TEST(view_to_simd_chunk_t, iterate)
{

}

TYPED_TEST(view_to_simd_chunk_t, with_padding)
{

}

TYPED_TEST(view_to_simd_chunk_t, without_padding)
{

}

TYPED_TEST(view_to_simd_chunk_t, nested)
{

}


TEST(to_simd_chunk_view, basic)
{
    using simd_t = simd_type_t<int32_t>;

    // for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
    //     std::ranges::copy(test::generate_sequence<dna4>(500, 10), std::back_inserter(data[i]));

    std::vector<dna4_vector> data1;
    for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
    {  // Use sequence generator.
                      // 0   1   2   3   4   5   6
        data1.push_back("ACGTAGCATCGACTGACGATCGACGACGC"_dna4);
        // data1.push_back("ACGTAGCATCGACTGACGATCGACGACG"_dna4);
        // data1.push_back("ACGTAGCATCGACTGACGATCGACGACG"_dna4);
        // data1.push_back("ACGTAGCATCGACTGACGATCGACGACG"_dna4);
    }

    // matrix:
// 0123|4567|8901|2345
//0ACGT|AGCA|TCGA|CTGA
//1ACGT|AGCA|TCGA|CTGA
//2ACGT|AGCA|TCGA|CTGA
//3ACGT|AGCA|TCGA|CTGA
//====================
//4CGAT|CGAC|GACG|----
//5CGAT|CGAC|GACG|----
//6CGAT|CGAC|GACG|----
//7CGAT|CGAC|GACG|----
//====================
//8----|----|----|----
//9----|----|----|----
//0----|----|----|----
//1----|----|----|----
//====================
//2----|----|----|----
//3----|----|----|----
//4----|----|----|----
//5----|----|----|----

// transposed

// 0123|4567|8901|2345
//0AAAA|CCCC|----|----
//1CCCC|GGGG|----|----
//2GGGG|AAAA|----|----
//3TTTT|TTTT|----|----
//====================
//4AAAA|CCCC|----|----
//5GGGG|GGGG|----|----
//6CCCC|AAAA|----|----
//7AAAA|CCCC|----|----
//====================
//8TTTT|GGGG|----|----
//9CCCC|AAAA|----|----
//0GGGG|CCCC|----|----
//1AAAA|GGGG|----|----
//====================
//2CCCC|----|----|----
//3TTTT|----|----|----
//4GGGG|----|----|----
//5AAAA|----|----|----

    // auto data = data1 | std::view::all;

    // detail::to_simd_chunk_view<decltype(data), simd_t> chunk_view{data, -1};

    auto chunk_view = data1 | view::to_simd_chunk<simd_t>;

    debug_stream << "Size: "<<  chunk_view.size() << "\n";

    size_t i = 0;
    for (auto & chunk : chunk_view)
    {
        debug_stream << "chunk " << i++ << ": ";
        for (auto & vec : chunk)
            debug_stream << vec << ", ";
        debug_stream << '\n';
    }
}
