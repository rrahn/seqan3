// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_debug_stream.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/range/view/to_simd_chunk.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(vectorised_alignment, smoke_test)
{
    using simd_t = simd_type_t<int16_t>;
    using simd_vec_t = std::vector<simd_t, aligned_allocator<simd_t, alignof(simd_t)>>;

    using cell_t = std::tuple<simd_t, simd_t, detail::ignore_t>;

    using result_t = alignment_result<typename align_result_selector<simd_vec_t, simd_vec_t, config_t>::type>;

    // These are the policies we need to define.
    using init_policy_t = deferred_crtp_base<affine_gap_init_policy>;
    using gap_policy_t  = deferred_crtp_base<affine_gap_policy, cell_t>;
    using opt_policy_t  = deferred_crtp_base<find_optimum_policy>;
    using mat_policy_t  = deferred_crtp_base<unbanded_score_dp_matrix_policy, aligned_allocator<cell_t, alignof(cell_t)>>;

    

    std::vector<dna4_vector> data1;
    for (size_t i = 0; i < 24; ++i)
    {  // Use sequence generator.
        data1.push_back("ACGTAGCATCGACTGACGATCGACGACT"_dna4);
    }

    std::vector<dna4_vector> data2;
    for (size_t i = 0; i < 24; ++i)
    {  // Use sequence generator.
        data2.push_back("GTAGCTAGCGACTAGGATAGATGCA"_dna4);
    }

    // This is our input, which we chunk on the simd_batch size.
    auto pairs = std::view::zip(data1, data2); //| view::chunk(simd_traits<simd_t>::length);

    // We need to transform the sequences.
    // Iterate over the batches
    // We would get a chunking range here.
    auto it = pairs.begin();
    auto s = pairs.end();
    size_t total_size = std::ranges::size(pairs); // we require the size here?
    size_t batch_offset = 0;

    // We iterate over the sequences.
    while (it != s)
    {
        // We need to take care of the last block.
        size_t batch_size = simd_traits<simd_t>::length;

        if (batch_size > total_size)
            batch_size -= batch_offset - total_size;

        debug_stream << "DEBUG: batch_size " << batch_size << " batch_offset " << batch_offset << "\n";
        auto pair_v = pairs | view::slice(batch_offset, batch_offset + batch_size);

        // Now we got the subrange which we have to transform to a simd vector.
        simd_vec_t simd_seq1;
        simd_vec_t simd_seq2;

        auto v1 = pair_v | view::get<0> | view::to_simd_chunk<simd_t>(0x8000) | std::view::join;
        auto v2 = pair_v | view::get<1> | view::to_simd_chunk<simd_t>(0x8000) | std::view::join;

        // Transform into simd batch vector.
        std::ranges::copy(v1, std::back_inserter(simd_seq1));
        std::ranges::copy(v2, std::back_inserter(simd_seq2));

        debug_stream << "DEBUG: " << simd_seq1[0] << " " << simd_seq2[0] << "\n";

        // next step: Call the alignment configurator and get alignment kernel back.



        // Move the iterator to the next chunk begin.
        std::ranges::advance(it, batch_size);
        batch_offset += batch_size;
    }

}
