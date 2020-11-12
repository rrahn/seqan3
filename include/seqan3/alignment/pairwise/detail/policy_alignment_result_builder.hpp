// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_result_builder.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment result builder.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration; must be a type specialisation of
 *                                   seqan3::configuration.
 *
 * \details
 *
 * Implements the interfaces to build the alignment result based on the previously selected output configurations.
 */
template <typename alignment_configuration_t>
#if !SEQAN3_WORKAROUND_GCC_93467
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
#endif // !SEQAN3_WORKAROUND_GCC_93467
class policy_alignment_result_builder
{
protected:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The alignment result type.
    using result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<result_type, empty_type>, "The alignment result type was not configured.");

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_result_builder() = default; //!< Defaulted.
    policy_alignment_result_builder(policy_alignment_result_builder const &) = default; //!< Defaulted.
    policy_alignment_result_builder(policy_alignment_result_builder &&) = default; //!< Defaulted.
    policy_alignment_result_builder & operator=(policy_alignment_result_builder const &) = default; //!< Defaulted.
    policy_alignment_result_builder & operator=(policy_alignment_result_builder &&) = default; //!< Defaulted.
    ~policy_alignment_result_builder() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration [not used in this context].
     */
    policy_alignment_result_builder(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {}
    //!\}


    /*!\brief Builds the seqan3::alignment_result based on the given alignment result type and then invokes the
     *        given callable with the result.
     *
     * \tparam logger_specs_t A template parameter pack over the template types of the
     *                        seqan3::detail::policy_alignment_algorithm_logger.
     * \tparam params_t A template parameter pack over the remaining parameters needed to construct the alignment
     *                  result.
     *
     * \param[in] logger The logger if applicable, otherwise the instance of the alignment algorithm.
     * \param[in] params The remaining parameters to construct the alignment result with.
     *
     * \details
     *
     * This is a wrapper function to promote the logger functionality of the alignment algorithm. The
     * seqan3::detail::policy_alignment_algorithm_logger is only available if the alignment was run in debug mode. In this case,
     * the seqan3::detail::policy_alignment_algorithm_logger is another base class of the alignment algorithm. When calling this
     * interface with the instance of the alignment algorithm, which is a derived class of the debug policy, then the
     * algorithm will be implicitly converted to the debug logger base class. Accordingly, we can access the public
     * members of the logger class inside of the result builder to store the logged alignment matrix.
     * This solves the problem of accessing the state of a sibling policy without weakening the interface of the
     * policies.
     *
     * If the debug mode was not enabled, another overload will be used, which does not promote the alignment matrix
     * to the debug logger class, i.e. calling it without debug information is still valid.
     */
    template <typename ...logger_specs_t, typename ...params_t>
    //!\cond
        requires traits_type::is_debug
    //!\endcond
    void make_result_and_invoke(policy_alignment_algorithm_logger<logger_specs_t...> & logger, params_t && ...params)
    {
        make_result_and_invoke_impl(logger, std::forward<params_t>(params)...);
    }

    //!\overload
    template <typename ...params_t>
    //!\cond
        requires (!traits_type::is_debug)
    //!\endcond
    void make_result_and_invoke(params_t && ...params)
    {
        make_result_and_invoke_impl(std::forward<params_t>(params)...);
    }

    /*!\brief Builds the seqan3::alignment_result based on the given alignment result type and then invokes the
     *        given callable with the result.
     *
     * \tparam logger_t The type the debug logger.
     * \tparam sequence_pair_t The type of the sequence pair.
     * \tparam id_t The type of the id.
     * \tparam score_t The type of the score.
     * \tparam matrix_coordinate_t The type of the matrix coordinate.
     * \tparam aligned_sequence_builder_t The type of the aligned sequences builder.
     * \tparam callback_t The type of the callback to invoke.
     *
     * \param[in] logger The debug logger (only available when run in debug mode).
     * \param[in] sequence_pair The indexed sequence pair.
     * \param[in] id The associated id.
     * \param[in] score The best alignment score.
     * \param[in] end_positions The matrix coordinate of the best alignment score.
     * \param[in] aligned_sequence_builder The aligned sequence builder.
     * \param[in] callback The callback to invoke with the generated result.
     *
     * \details
     *
     * Generates a seqan3::alignment_result object with the results computed during the alignment. Depending on the
     * \ref seqan3_align_cfg_output_configurations "seqan3::align_cfg::output_*" configuration only the requested values
     * are stored. In some cases some additional work is done to generate the requested result. For example computing
     * the associated alignment from the traceback matrix.
     *
     * The first parameter is the logger and is only available if the alignment was run in debug mode. It stores the
     * debug score and, if applicable, the debug trace matrix, which are then stored in the created alignment result.
     */
    template <typename logger_t,
              typename sequence_pair_t,
              typename index_t,
              typename score_t,
              typename matrix_coordinate_t,
              typename aligned_sequence_builder_t,
              typename callback_t>
    //!\cond
        requires std::invocable<callback_t, result_type>
    //!\endcond
    void make_result_and_invoke_impl([[maybe_unused]] logger_t && logger,
                                     [[maybe_unused]] sequence_pair_t && sequence_pair,
                                     [[maybe_unused]] index_t && id,
                                     [[maybe_unused]] score_t score,
                                     [[maybe_unused]] matrix_coordinate_t end_positions,
                                     [[maybe_unused]] aligned_sequence_builder_t const & aligned_sequence_builder,
                                     callback_t && callback)
    {
        using std::get;
        using invalid_t = std::nullopt_t *;

        result_type result{};

        if constexpr (traits_type::output_sequence1_id)
            result.data.sequence1_id = id;

        if constexpr (traits_type::output_sequence2_id)
            result.data.sequence2_id = id;

        if constexpr (traits_type::compute_score)
        {
            static_assert(!std::same_as<decltype(result.data.score), invalid_t>,
                          "Invalid configuration. Expected result with score!");
            result.data.score = std::move(score);
        }

        if constexpr (traits_type::compute_end_positions)
        {
            static_assert(!std::same_as<decltype(result.data.end_positions), invalid_t>,
                          "Invalid configuration. Expected result with end positions!");

            result.data.end_positions.first = end_positions.col;
            result.data.end_positions.second = end_positions.row;
        }

        if constexpr (traits_type::requires_trace_information)
        {
            auto aligned_sequence_result = aligned_sequence_builder(get<0>(sequence_pair),
                                                                    get<1>(sequence_pair),
                                                                    end_positions);

            if constexpr (traits_type::compute_begin_positions)
            {
                result.data.begin_positions.first = aligned_sequence_result.first_sequence_slice_positions.first;
                result.data.begin_positions.second = aligned_sequence_result.second_sequence_slice_positions.first;
            }

            if constexpr (traits_type::compute_sequence_alignment)
                result.data.alignment = std::move(aligned_sequence_result.alignment);
        }

        // In case we run in debug mode, we need to store the score and possibly trace matrix.
        if constexpr (traits_type::is_debug)
        {
            // static constexpr that logger is not a ignore type.
            using std::swap;
            swap(result.data.score_debug_matrix, logger.debug_score_matrix);

            if (traits_type::compute_sequence_alignment)
                swap(result.data.trace_debug_matrix, logger.debug_trace_matrix);
        }

        callback(std::move(result));
    }

};
} // namespace seqan3::detail
