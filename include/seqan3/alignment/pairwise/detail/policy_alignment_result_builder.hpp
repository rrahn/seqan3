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

#include <seqan3/std/concepts>

#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment result builder.
 * \ingroup pairwise_alignment
 *
 * \tparam result_t The pairwise alignment result type; must model seqan3::detail::template_specialisation_of
 *                  seqan3::alignment_result.
 *
 * \details
 *
 * Implements the interfaces to build the alignment result based on the previously selected output configurations.
 */
template <typename result_t>
//!\cond
#if !SEQAN3_WORKAROUND_GCC_93467
    requires is_type_specialisation_of_v<result_t, alignment_result>
#endif // !SEQAN3_WORKAROUND_GCC_93467
//!\endcond
class policy_alignment_result_builder
{
protected:
    // ----------------------------------------------------------------------------
    // type defintions
    // ----------------------------------------------------------------------------

    //!\brief The alignment result type.
    using result_type = result_t;

    static_assert(!std::same_as<result_type, empty_type>, "The alignment result type was not configured.");

    //!\brief The internal data type storing the alignment information.
    using result_data_type = typename alignment_result_value_type_accessor<result_type>::type;

    // ----------------------------------------------------------------------------
    // static member
    // ----------------------------------------------------------------------------

    //!\brief Flag indicating whether the trace matrix was computed.
    static constexpr bool requires_trace_information = result_data_type::has_begin_positions ||
                                                       result_data_type::has_alignment;
    //!\brief Flag indicating whether the debug mode was enabled.
    static constexpr bool is_debug = result_data_type::has_debug_score_matrix ||
                                     result_data_type::has_debug_trace_matrix;

    // ----------------------------------------------------------------------------
    // member function
    // ----------------------------------------------------------------------------

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
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
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
        requires is_debug
    //!\endcond
    void make_result_and_invoke(policy_alignment_algorithm_logger<logger_specs_t...> & logger, params_t && ...params)
    {
        make_result_and_invoke_impl(logger, std::forward<params_t>(params)...);
    }

    //!\overload
    template <typename ...params_t>
    //!\cond
        requires (!is_debug)
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
     * \tparam trace_path_from_t The type of the aligned sequences builder.
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
        result_type result{};
        result_data_type & result_data = alignment_result_attorney::data(result);

        if constexpr (result_data_type::has_sequence1_id)
            result_data.sequence1_id = id;

        if constexpr (result_data_type::has_sequence2_id)
            result_data.sequence2_id = id;

        if constexpr (result_data_type::has_score)
            result_data.score = std::move(score);

        if constexpr (result_data_type::has_end_positions)
        {
            result_data.end_positions.first = end_positions.col;
            result_data.end_positions.second = end_positions.row;
        }

        if constexpr (requires_trace_information)
        {
            using std::get;

            auto aligned_sequence_result = aligned_sequence_builder(get<0>(sequence_pair),
                                                                    get<1>(sequence_pair),
                                                                    end_positions);

            if constexpr (result_data_type::has_begin_positions)
            {
                result_data.begin_positions.first = aligned_sequence_result.first_sequence_slice_positions.first;
                result_data.begin_positions.second = aligned_sequence_result.second_sequence_slice_positions.first;
            }

            if constexpr (result_data_type::has_alignment)
                result_data.alignment = std::move(aligned_sequence_result.alignment);
        }

        // In case we run in debug mode, we need to store the score and possibly trace matrix.
        if constexpr (is_debug)
        {
            // static constexpr that logger is not a ignore type.
            using std::swap;
            swap(result_data.score_debug_matrix, logger.debug_score_matrix);

            if (result_data_type::has_alignment)
                swap(result_data.trace_debug_matrix, logger.debug_trace_matrix);
        }

        callback(std::move(result));
    }
};
} // namespace seqan3::detail
