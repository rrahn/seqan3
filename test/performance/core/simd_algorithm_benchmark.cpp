// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <array>
#include <cstdlib>

#include <benchmark/benchmark.h>

#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

inline auto make_matrix()
{
    using simd_t = simd_type_t<int8_t>;

    std::array<simd_t, simd_traits<simd_t>::length> matrix;
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            matrix[i][j] = std::rand() % 10;

    return matrix;
}

template <typename simd_t>
inline auto reduce(simd_t const & vec)
{
    size_t sum = 0;
    for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        sum += vec[i];

    return sum;
}

// ----------------------------------------------------------------------------
// Benchhmark transpose
// ----------------------------------------------------------------------------

template <typename simd_t>
inline void generic_transpose(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    std::array<simd_t, simd_traits<simd_t>::length> tmp;
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            tmp[j][i] = matrix[i][j];

    std::swap(tmp, matrix);
}

template <bool use_intrinsics>
static void transpose(benchmark::State& state)
{
    size_t sum = 0;

    for (auto _ : state)
    {
        for (size_t i = 0; i < 100; ++i)
        {
            state.PauseTiming();
            auto matrix = make_matrix();
            state.ResumeTiming();

            if constexpr (use_intrinsics)
                simd::transpose(matrix);
            else
                generic_transpose(matrix);

            state.PauseTiming();
            sum += reduce(matrix[std::rand() % matrix.size()]);
            state.ResumeTiming();
        }
    }

    state.counters["checksum"] = sum;
}

static void transpose_intrinsics(benchmark::State& state)
{
    transpose<true>(state);
}

static void transpose_autovec(benchmark::State& state)
{
    transpose<false>(state);
}

BENCHMARK(transpose_intrinsics);
BENCHMARK(transpose_autovec);

template <bool use_simd_upcast, typename source_t, typename target_t>
static void do_upcast(benchmark::State& state)
{
    source_t src = simd::iota<source_t>(std::rand() % 100);
    target_t target{};
    size_t sum = 0;

    for (auto _ : state)
    {
        for (size_t i = 0; i < 1'000; ++i)
        {
            if constexpr (use_simd_upcast)
                target = simd::upcast<target_t>(src);
            else
                target = simd::upcast_auto<target_t>(src);

            state.PauseTiming();
            sum += reduce(target);
            state.ResumeTiming();
        }
    }

    state.counters["checksum"] = sum;
}

// ----------------------------------------------------------------------------
// Benchhmark upcast
// ----------------------------------------------------------------------------

template <typename source_t, typename target_t>
static void upcast_intrinsics(benchmark::State& state)
{
    do_upcast<true, source_t, target_t>(state);
}

template <typename source_t, typename target_t>
static void upcast_autovec(benchmark::State& state)
{
    do_upcast<false, source_t, target_t>(state);
}

BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int8_t>, simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int8_t>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int8_t>, simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int16_t>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int16_t>, simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast_intrinsics, simd_type_t<int32_t>, simd_type_t<int64_t>);

BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int8_t>, simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int8_t>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int8_t>, simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int16_t>, simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int16_t>, simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast_autovec, simd_type_t<int32_t>, simd_type_t<int64_t>);

BENCHMARK_MAIN();
