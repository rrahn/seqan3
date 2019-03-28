// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/simd_scoring_scheme_simple.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;
using namespace seqan3::detail;

TEST(simd_scoring_scheme_wrapper, construction)
{
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, dna4>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<scheme_t>);
    EXPECT_TRUE((std::is_constructible_v<scheme_t, nucleotide_scoring_scheme<int16_t>>));
    EXPECT_TRUE(std::Semiregular<scheme_t>);
}

TEST(simd_scoring_scheme_wrapper, construct_from_scoring_scheme)
{
    using simd_t = typename simd::simd_type<int8_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, dna4>;

    { // no throw
        scheme_t scheme{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}};

        simd_t s1 = simd::fill<simd_t>(2);
        simd_t s2 = simd::fill<simd_t>(2);
        SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), simd::fill<simd_t>(4));

        s2 = simd::fill<simd_t>(1);
        SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), simd::fill<simd_t>(-5));
    }

    { // throw if values would overflow
        EXPECT_THROW((scheme_t{nucleotide_scoring_scheme<int16_t>{match_score{128}, mismatch_score{-5}}}),
                     std::invalid_argument);
        EXPECT_THROW((scheme_t{nucleotide_scoring_scheme<int16_t>{match_score{4}, mismatch_score{-129}}}),
                     std::invalid_argument);
    }
}

TEST(simd_scoring_scheme_wrapper, score_global)
{
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, dna4>;

    scheme_t scheme{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}};

    simd_t s1 = simd::fill<simd_t>(2);
    simd_t s2 = simd::fill<simd_t>(2);

    // all match
    SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), simd::fill<simd_t>(4));

    // all mismatch
    s2 = simd::fill<simd_t>(3);
    SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), simd::fill<simd_t>(-5));

    // first one is padded so it must be a match.
    s2[0] = 1 << ((sizeof(typename simd_traits<simd_t>::scalar_type) << 3) - 1);
    simd_t res = simd::fill<simd_t>(-5);
    res[0] = 4;
    SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), res);

    // first one in other sequence is padded so it must be still a match.
    s1[0] = 1 << ((sizeof(typename simd_traits<simd_t>::scalar_type) << 3) - 1);
    SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), res);

    // Only first one in other sequence is padded so it must be still a match.
    s2[0] = 3;
    SIMD_EQ(scheme.score<global_alignment_type>(s1, s2), res);
}

TEST(simd_scoring_scheme_wrapper, score_local)
{
    // In local alignment we always want to mismatch.
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, dna4>;

    scheme_t scheme{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}};

    simd_t s1 = simd::fill<simd_t>(2);
    simd_t s2 = simd::fill<simd_t>(2);

    // all match
    SIMD_EQ(scheme.score<local_alignment_type>(s1, s2), simd::fill<simd_t>(4));

    // all mismatch
    s2 = simd::fill<simd_t>(3);
    SIMD_EQ(scheme.score<local_alignment_type>(s1, s2), simd::fill<simd_t>(-5));

    // first one is padded so it must be a mismatch.
    s2 = simd::fill<simd_t>(2);
    s2[0] = 1 << ((sizeof(typename simd_traits<simd_t>::scalar_type) << 3) - 1);
    simd_t res = simd::fill<simd_t>(4);
    res[0] = -5;
    SIMD_EQ(scheme.score<local_alignment_type>(s1, s2), res);

    // first one in other sequence is padded as well so it must be still a mismatch.
    s1[0] = 1 << ((sizeof(typename simd_traits<simd_t>::scalar_type) << 3) - 2);
    SIMD_EQ(scheme.score<local_alignment_type>(s1, s2), res);

    // Only first one in other sequence is padded so it must be still a mismatch.
    s2[0] = 3;
    SIMD_EQ(scheme.score<local_alignment_type>(s1, s2), res);
}
