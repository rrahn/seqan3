// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/detail/select_alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/detail/select_alignment_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/detail/select_recursion_policy.hpp>
#include <seqan3/alignment/pairwise/detail/select_result_builder_policy.hpp>
#include <seqan3/alignment/pairwise/detail/select_logger_policy.hpp>
#include <seqan3/alignment/pairwise/detail/select_optimum_tracker_policy.hpp>
#include <seqan3/alignment/pairwise/detail/select_scoring_scheme_policy.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>

namespace fixture = seqan3::test::alignment::fixture;

/*!\brief The alignment base fixture.
 *
 * \tparam _fixture A non-type pointer template argument that points to the test fixture for this test.
 *
 * \details
 *
 * The non-type template argument allows to output a diagnostic message that contains the fully qualified name of the
 * used test fixture when a test case fails. This helps readers of the test output to identify the exact location
 * of the data used in the failed test.
 */
template <auto _fixture>
struct pairwise_alignment_fixture : public ::testing::Test
{
    //!\brief Return the fixture as concrete instance (to avoid pointer-semantics).
    auto fixture() -> decltype(fixture::alignment_fixture_collection{*_fixture}) const &
    {
        return *_fixture;
    }
};

/*!\brief Type trait to generate the alignment for the given configuration and the alphabet type.
 *
 * \tparam config_t The alignment configuration type.
 * \tparam alphabet_t The sequence alphabet type.
 */
template <typename config_t, typename alphabet_t>
struct alignment_algorithm_for
{
private:
    using traits_t = seqan3::detail::alignment_configuration_traits<config_t>;
    using optimum_tracker_policy_t = seqan3::detail::select_optimum_policy_t<traits_t>;
    using gap_cost_policy_t = seqan3::detail::select_recursion_policy_t<traits_t>;
    using result_builder_policy_t = seqan3::detail::select_result_builder_policy_t<traits_t>;
    using scoring_scheme_policy_t = seqan3::detail::select_scoring_scheme_policy_t<traits_t, alphabet_t>;
    using alignment_matrix_policy_t = seqan3::detail::select_alignment_matrix_policy_t<traits_t>;
    using alignment_logger_policy_t = seqan3::detail::select_logger_policy_t<traits_t>;

public:
    using type = seqan3::detail::select_alignment_algorithm_t<traits_t,
                                                              gap_cost_policy_t,
                                                              scoring_scheme_policy_t,
                                                              optimum_tracker_policy_t,
                                                              alignment_logger_policy_t,
                                                              alignment_matrix_policy_t,
                                                              result_builder_policy_t>;
};

template <typename fixture_t>
class simd_pairwise_alignment_test : public fixture_t
{
    // Define sequence types.
    using alphabet_t = seqan3::dna4;
    using sequence_t = std::vector<alphabet_t>;
    using sequence_pair_t = std::pair<sequence_t, sequence_t>;
    using indexed_sequence_pair_t = std::pair<sequence_pair_t, size_t>;
    using indexed_sequence_pairs_t = std::vector<indexed_sequence_pair_t>;
    using chunked_sequence_pairs_t = std::vector<indexed_sequence_pairs_t>;

    chunked_sequence_pairs_t chunked_sequences{};

    //!\brief Initialise the sequences for the current test case run.
    void SetUp() override
    {
        auto const & fixture = this->fixture();

        using configuration_t = decltype(fixture.config);

        static_assert(configuration_t::template exists<seqan3::align_cfg::score_type>(),
                      "The test template requires the score type. Please provide seqan3::align_cfg::score_type "
                      "in your test type definitions!");

        using score_config_t =
            std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::score_type>(fixture.config))>;
        using score_t = typename score_config_t::type;

        // Buid the indexed sequences for the test case.
        indexed_sequence_pairs_t indexed_sequences{};
        chunked_sequences.clear();

        size_t const sequence_count = fixture.collection.size();
        constexpr size_t batch_size = seqan3::simd_traits<seqan3::simd_type_t<score_t>>::length;
        for (size_t seq_idx = 0; seq_idx < sequence_count;)
        {
            indexed_sequences.clear();
            for (size_t batch_count = 0; seq_idx < sequence_count && batch_count < batch_size; ++seq_idx, ++batch_count)
            {
                indexed_sequences.push_back(
                        indexed_sequence_pair_t{sequence_pair_t{fixture.collection[seq_idx].sequence1,
                                                                fixture.collection[seq_idx].sequence2},
                                                seq_idx});
            }
            chunked_sequences.push_back(indexed_sequences);
        }
    }

    /*!\brief Configures the alignment result based on the given configuration object.
     *
     * \tparam config_t The configuration type.
     *
     * \param[in] config The configuration to add the alignment result type.
     *
     * \returns A new configuration with added alignment result type.
     */
    template <typename config_t>
    auto configure_alignment_result(config_t const & config)
    {
        using alignment_result_value_t = typename seqan3::detail::align_result_selector<sequence_t &, sequence_t &, config_t>::type;
        using alignment_result_t = seqan3::alignment_result<alignment_result_value_t>;

        return config | seqan3::align_cfg::detail::result_type<alignment_result_t>{};
    }

public:
    /*!\brief Computes the alignment for the given alignment configuration.
     *
     * \tparam config_t The configuration type.
     * \tparam callback_t The type of the callback invoked for every computed alignment result.
     *
     * \param[in] config The alignment configuration.
     * \param[in] callback The test callback invoked on every alignment result.
     */
    template <typename config_t, typename callback_t>
    void compute_alignment(config_t const & config, callback_t && callback)
    {
        auto complete_config = configure_alignment_result(config);
        using algorithm_t = typename alignment_algorithm_for<decltype(complete_config), alphabet_t>::type;

        // We want to further configure the alignments here.
        algorithm_t alignment{complete_config};
        for (auto & chunk : chunked_sequences)
            alignment(chunk, callback);
    }
};

TYPED_TEST_SUITE_P(simd_pairwise_alignment_test);

TYPED_TEST_P(simd_pairwise_alignment_test, score)
{
    auto expected_scores = this->fixture().get_scores();
    size_t expected_index = 0;
    auto config = this->fixture().config | seqan3::align_cfg::vectorised{} | seqan3::align_cfg::output_score{};
    this->compute_alignment(config, [&](auto res)
    {
        EXPECT_EQ(res.score(), expected_scores[expected_index]) << "at index: " << expected_index;
        ++expected_index;
    });
}

TYPED_TEST_P(simd_pairwise_alignment_test, end_position)
{
    auto expected_end_positions = this->fixture().get_end_positions();
    size_t expected_index = 0;
    auto config = this->fixture().config | seqan3::align_cfg::vectorised{} | seqan3::align_cfg::output_end_position{};

    this->compute_alignment(config, [&](auto res)
    {
        EXPECT_EQ(res.sequence1_end_position(), expected_end_positions[expected_index].first)
            << "at index: " << expected_index;
        EXPECT_EQ(res.sequence2_end_position(), expected_end_positions[expected_index].second)
            << "at index: " << expected_index;
        ++expected_index;
    });
}

TYPED_TEST_P(simd_pairwise_alignment_test, begin_position)
{
    auto expected_begin_positions = this->fixture().get_begin_positions();
    size_t expected_index = 0;
    auto config = this->fixture().config | seqan3::align_cfg::vectorised{} | seqan3::align_cfg::output_begin_position{};

    this->compute_alignment(config, [&](auto res)
    {
        EXPECT_EQ(res.sequence1_begin_position(), expected_begin_positions[expected_index].first)
            << "at index: " << expected_index;
        EXPECT_EQ(res.sequence2_begin_position(), expected_begin_positions[expected_index].second)
            << "at index: " << expected_index;
        ++expected_index;
    });
}

TYPED_TEST_P(simd_pairwise_alignment_test, alignment)
{
    auto expected_aligned1 = this->fixture().get_aligned_sequences1();
    auto expected_aligned2 = this->fixture().get_aligned_sequences2();
    size_t expected_index = 0;
    auto config = this->fixture().config | seqan3::align_cfg::vectorised{} | seqan3::align_cfg::output_alignment{};

    this->compute_alignment(config, [&](auto res)
    {
        using std::get;
        EXPECT_RANGE_EQ(get<0>(res.alignment()) | seqan3::views::to_char, expected_aligned1[expected_index]);
        EXPECT_RANGE_EQ(get<1>(res.alignment()) | seqan3::views::to_char, expected_aligned2[expected_index]);
        ++expected_index;
    });
}

REGISTER_TYPED_TEST_SUITE_P(simd_pairwise_alignment_test, score, end_position, begin_position, alignment);
