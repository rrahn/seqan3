// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "../fixture/global_affine_banded.hpp"
#include "simd_alignment_test_template.hpp"

namespace test_fixture = seqan3::test::alignment::fixture;

static auto dna4_all_same = []()
{
    auto base_fixture = test_fixture::global::affine::banded::dna4_01;
    using fixture_t = decltype(base_fixture);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 100; ++i)
        data.push_back(base_fixture);

    return test_fixture::alignment_fixture_collection{base_fixture.config | seqan3::align_cfg::score_type<int32_t>{},
                                                      data};
}();

// static auto dna4_different_length = []()
// {
//     auto base_fixture_01 = test_fixture::global::affine::banded::dna4_match_4_mismatch_5_gap_1_open_10_part_01;
//     auto base_fixture_02 = test_fixture::global::affine::banded::dna4_match_4_mismatch_5_gap_1_open_10_part_02;
//     auto base_fixture_03 = test_fixture::global::affine::banded::dna4_match_4_mismatch_5_gap_1_open_10_part_03;
//     auto base_fixture_04 = test_fixture::global::affine::banded::dna4_match_4_mismatch_5_gap_1_open_10_part_04;

//     using fixture_t = decltype(base_fixture_01);

//     std::vector<fixture_t> data{};
//     for (size_t i = 0; i < 25; ++i)
//     {
//         data.push_back(base_fixture_01);
//         data.push_back(base_fixture_02);
//         data.push_back(base_fixture_03);
//         data.push_back(base_fixture_04);
//     }

//     return test_fixture::alignment_fixture_collection{base_fixture_01.config | seqan3::align_cfg::score_type<int32_t>{},
//                                                       data};
// }();

using test_types = ::testing::Types<pairwise_alignment_fixture<&dna4_all_same>>;
INSTANTIATE_TYPED_TEST_SUITE_P(simd_global_affine_banded_dna_test, simd_pairwise_alignment_test, test_types, );
