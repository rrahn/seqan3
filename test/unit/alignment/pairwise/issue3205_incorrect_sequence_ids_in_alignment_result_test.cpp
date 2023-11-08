// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

using seqan3::operator""_dna4;

TEST(issue2305, incorrect_sequence_ids)
{
    std::vector vec{"ACGTGACTGACT"_dna4, "ACGAAGACCGAT"_dna4, "ACGTGACTGACT"_dna4, "AGGTACGAGCGACACT"_dna4};

    auto config = seqan3::align_cfg::method_global{} |
                  seqan3::align_cfg::edit_scheme |
                  seqan3::align_cfg::min_score{-7} |
                  seqan3::align_cfg::output_score{} |
                  seqan3::align_cfg::output_sequence1_id{} |
                  seqan3::align_cfg::output_sequence2_id{};

    auto index_pairs = seqan3::views::pairwise_combine(std::views::iota(0ul, vec.size()));
    auto alignment_results = seqan3::align_pairwise(seqan3::views::pairwise_combine(vec), config);
    auto filter_v = std::views::filter( [](auto && res) { return res.score() >= -6 ;  });

    for (auto const & result : alignment_results | filter_v)
        seqan3::debug_stream << index_pairs[result.sequence1_id()] << '\n';
}
