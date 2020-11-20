// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/std/ranges>
#include <tuple>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
    #include <seqan/align_parallel.h>
#endif

// Globally defined constants to ensure same test data.
inline constexpr size_t sequence_length = 150;
#ifndef NDEBUG
inline constexpr size_t set_size        = 128;
#else
inline constexpr size_t set_size        = 1024;
#endif // NDEBUG

// We don't know if the system supports hyper-threading so we use only half the threads so that the
// simd benchmark is likely to run on physical cores only.
uint32_t get_number_of_threads()
{
    uint32_t thread_count = std::thread::hardware_concurrency();
    return (thread_count == 1) ? thread_count : thread_count >> 1;
}

namespace sd = seqan3::detail;

template <typename config_t, typename alphabet_t>
struct alignment_algorithm_for
{
private:
    using traits_t = sd::alignment_configuration_traits<config_t>;
    using optimum_tracker_policy_t = sd::select_optimum_policy_t<traits_t>;
    using gap_cost_policy_t = sd::select_recursion_policy_t<traits_t>;
    using result_builder_policy_t = sd::select_result_builder_policy_t<traits_t>;
    using scoring_scheme_policy_t = sd::select_scoring_scheme_policy_t<traits_t, alphabet_t>;
    using alignment_matrix_policy_t = sd::select_alignment_matrix_policy_t<traits_t>;
    using alignment_logger_policy_t = sd::select_logger_policy_t<traits_t>;
public:
    using type = sd::select_alignment_algorithm_t<traits_t,
                                                  gap_cost_policy_t,
                                                  scoring_scheme_policy_t,
                                                  optimum_tracker_policy_t,
                                                  alignment_logger_policy_t,
                                                  alignment_matrix_policy_t,
                                                  result_builder_policy_t>;
};

template <typename config_t, typename data_t>
auto configure_alignment_result(config_t const & config, data_t const &)
{
    using sequence_pair_t = std::ranges::range_value_t<data_t>;
    using sequence_t = std::tuple_element_t<0, sequence_pair_t>;
    using alignment_result_value_t = typename sd::align_result_selector<sequence_t &, sequence_t &, config_t>::type;
    using alignment_result_t = seqan3::alignment_result<alignment_result_value_t>;

    return config | seqan3::align_cfg::detail::result_type<alignment_result_t>{};
}


// ----------------------------------------------------------------------------
//  seqan3 pairwise alignment
// ----------------------------------------------------------------------------

template <typename alphabet_t, typename ...align_configs_t>
void seqan3_affine_accelerated(benchmark::State & state, alphabet_t, align_configs_t && ...configs)
{
    size_t sequence_length_variance = state.range(0);
    auto data = seqan3::test::generate_sequence_pairs<alphabet_t>(sequence_length,
                                                                  set_size,
                                                                  sequence_length_variance);

    int64_t total = 0;
    auto accelerate_config = configure_alignment_result((configs | ...), data);

    using algo_t = typename alignment_algorithm_for<decltype(accelerate_config), alphabet_t>::type;
    algo_t algo{accelerate_config};

    using score_t = typename std::remove_reference_t<decltype(accelerate_config.get_or(seqan3::align_cfg::score_type<int32_t>{}))>::type;

    constexpr size_t chunk_size = seqan3::simd_traits<seqan3::simd_type_t<score_t>>::length;
    auto chunked_data = seqan3::views::zip(data, std::views::iota(0)) | seqan3::views::chunk(chunk_size);

    for (auto _ : state)
    {
        for (auto && chunk : chunked_data)
            algo(chunk, [&] (auto res) { total += res.score(); });
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(data, accelerate_config);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#ifdef SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
//  seqan2 pairwise alignment
// ----------------------------------------------------------------------------

template <typename alphabet_t, typename ...args_t>
void seqan2_affine_accelerated(benchmark::State & state, alphabet_t, args_t && ...args)
{
    std::tuple captured_args{args...};

    size_t sequence_length_variance = state.range(0);
    auto [vec1, vec2] = seqan3::test::generate_sequence_pairs_seqan2<alphabet_t>(sequence_length,
                                                                                 set_size,
                                                                                 sequence_length_variance);

    auto scoring_scheme = std::get<0>(captured_args);
    auto exec = std::get<1>(captured_args);
    setNumThreads(exec, std::get<2>(captured_args));

    // Possibly enable banded alignment.
    auto seqan3_align_cfg = std::get<3>(captured_args);
    static constexpr bool execute_with_band = seqan3_align_cfg.template exists<seqan3::align_cfg::band_fixed_size>();
    [[maybe_unused]] int lower_diagonal{};
    [[maybe_unused]] int upper_diagonal{};

    if constexpr (execute_with_band)
    {
        using std::get;
        auto band_cfg = get<seqan3::align_cfg::band_fixed_size>(seqan3_align_cfg);
        lower_diagonal = band_cfg.lower_diagonal;
        upper_diagonal = band_cfg.upper_diagonal;
    }

    int64_t total = 0;
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::String<int> res;
        if constexpr (execute_with_band)
            res = seqan::globalAlignmentScore(exec, vec1, vec2, scoring_scheme, lower_diagonal, upper_diagonal);
        else
            res = seqan::globalAlignmentScore(exec, vec1, vec2, scoring_scheme);

        total = std::accumulate(seqan::begin(res), seqan::end(res), total);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), seqan3_align_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#endif // SEQAN3_HAS_SEQAN2
