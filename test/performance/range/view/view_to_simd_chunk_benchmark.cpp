// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <iterator>
#include <list>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/range/view/to_simd_chunk.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/performance/sequence_generator.hpp>

using namespace seqan3;

// ============================================================================
//  sequential_read
// ============================================================================

template <typename container_t, typename simd_t>
void to_simd(benchmark::State& state)
{
    // Preparing the sequences
    std::vector<container_t> data;
    data.resize(simd_traits<simd_t>::length);

    for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        std::ranges::copy(test::generate_sequence<dna4>(500, 10), std::back_inserter(data[i]));

    auto v = data | view::to_simd_chunk<simd_t>;
    size_t value = 0;
    for (auto _ : state)
    {
        size_t idx = std::rand() % simd_traits<simd_t>::length;
        for (auto & chunk : v)
            for (simd_t const & vec : chunk)
                value += vec[idx];
    }

    state.counters["value"] = value;
}

// runs with ContiguousRange
BENCHMARK_TEMPLATE(to_simd, std::vector<dna4>, simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<dna4>, simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<dna4>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<dna4>, simd_type_t<int64_t>);

BENCHMARK_TEMPLATE(to_simd, std::deque<dna4>, simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<dna4>, simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<dna4>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<dna4>, simd_type_t<int64_t>);

// runs without ContiguousRange
BENCHMARK_TEMPLATE(to_simd, std::list<dna4>, simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<dna4>, simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<dna4>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<dna4>, simd_type_t<int64_t>);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
