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

#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
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
    //!\brief The debug score type.
    using debug_score_type = std::optional<typename traits_type::score_type>;
    //!\brief The debug trace type.
    using debug_trace_type = std::optional<typename traits_type::trace_type>;
    //!\brief The type of the debug score matrix.
    using debug_score_matrix_t =
        std::conditional_t<traits_type::is_debug,
                           two_dimensional_matrix<debug_score_type,
                                                  std::allocator<debug_score_type>,
                                                  matrix_major_order::column>,
                           empty_type>;
    //!\brief The type of the debug trace matrix.
    using debug_trace_matrix_t =
        std::conditional_t<traits_type::is_debug && traits_type::compute_sequence_alignment,
                           two_dimensional_matrix<debug_trace_type,
                                                  std::allocator<debug_trace_type>,
                                                  matrix_major_order::column>,
                           empty_type>;

    static_assert(!std::same_as<result_type, empty_type>, "The alignment result type was not configured.");

    //!\brief The debug score matrix.
    debug_score_matrix_t debug_score_matrix{};
    //!\brief The debug trace matrix.
    debug_trace_matrix_t debug_trace_matrix{};

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
     * \tparam sequence_pair_t The type of the sequence pair.
     * \tparam id_t The type of the id.
     * \tparam score_t The type of the score.
     * \tparam matrix_coordinate_t The type of the matrix coordinate.
     * \tparam alignment_matrix_t The type of the alignment matrix.
     * \tparam callback_t The type of the callback to invoke.
     *
     * \param[in] sequence_pair The indexed sequence pair.
     * \param[in] id The associated id.
     * \param[in] score The best alignment score.
     * \param[in] end_positions The matrix coordinate of the best alignment score.
     * \param[in] alignment_matrix The alignment matrix to obtain the trace back from.
     * \param[in] callback The callback to invoke with the generated result.
     *
     * \details
     *
     * Generates a seqan3::alignment_result object with the results computed during the alignment. Depending on the
     * seqan3::align_cfg::result configuration only the requested values are stored. In some cases some additional
     * work is done to generate the requested result. For example computing the associated alignment from the traceback
     * matrix.
     */
    template <typename sequence_pair_t,
              typename index_t,
              typename score_t,
              typename matrix_coordinate_t,
              typename alignment_matrix_t,
              typename callback_t>
    //!\cond
        requires std::invocable<callback_t, result_type>
    //!\endcond
    void make_result_and_invoke([[maybe_unused]] sequence_pair_t && sequence_pair,
                                [[maybe_unused]] index_t && id,
                                [[maybe_unused]] score_t score,
                                [[maybe_unused]] matrix_coordinate_t end_positions,
                                [[maybe_unused]] alignment_matrix_t const & alignment_matrix,
                                callback_t && callback)
    {
        using std::get;
        using invalid_t = std::nullopt_t *;

        result_type result{};

        static_assert(!std::same_as<decltype(result.data.id), invalid_t>,
                      "Invalid configuration. Expected result with id!");
        result.data.id = std::move(id);

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
            aligned_sequence_builder builder{get<0>(sequence_pair), get<1>(sequence_pair)};
            auto aligned_sequence_result = builder(alignment_matrix.trace_path(end_positions));

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
            result.data.score_debug_matrix = debug_score_matrix;

            if (traits_type::compute_sequence_alignment)
                result.data.trace_debug_matrix = debug_trace_matrix;
        }

        callback(std::move(result));
    }

    /*!\brief Initialises the local debug matrices.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     *
     * \details
     *
     * Resizes the debug score and if requested the trace matrix to the given matrix dimensions.
     *
     * ### Exception
     *
     * Might throw std::bad_alloc if the requested matrix size exceeds the available memory.
     */
    void initialise_debug_matrices(size_t const sequence1_size, size_t const sequence2_size)
    {
        assert(sequence1_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        assert(sequence2_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));

        size_t const column_count = sequence1_size + 1;
        size_t const row_count = sequence2_size + 1;

        debug_score_matrix.resize(number_rows{row_count}, number_cols{column_count});

        if constexpr (traits_type::compute_sequence_alignment)
            debug_trace_matrix.resize(number_rows{row_count}, number_cols{column_count});
    }

    /*!\brief Log the current alignment column.
     * \tparam coordinate_column_t The type of the coordinate matrix column; must model std::ranges::input_range and
     *                             seqan3::detail::matrix_offset is std::constructible_from the range reference type.
     * \tparam alignment_column_t The type of the alignment matrix column; must model std::ranges::input_range and
     *                            the range reference type must model seqan3::detail::affine_cell_proxy_instance.
     *
     * param[in] coordinate_column The current column over the coordinate matrix.
     * param[in] alignment_column The current column over the alignment matrix.
     *
     * \details
     *
     * Logs the current alignment column in the locally stored debug matrices for the score and trace matrix.
     * The coordinate column is used to store the column at the correct offset. This is needed when logging
     * the banded matrix, where the column offset can be different. The trace matrix is only logged if the
     * sequence alignment shall be computed.
     */
    template <std::ranges::input_range coordinate_column_t, std::ranges::input_range alignment_column_t>
    //!\cond
        requires (std::constructible_from<matrix_offset, std::ranges::range_reference_t<coordinate_column_t>> &&
                  affine_cell_proxy_instance<std::ranges::range_reference_t<alignment_column_t>>)
    //!\endcond
    void log_alignment_matrix_column(coordinate_column_t && coordinate_column,
                                     alignment_column_t && alignment_column)
    {
        matrix_offset column_coordinate_begin{*coordinate_column.begin()};

        assert(static_cast<size_t>(column_coordinate_begin.col) < debug_score_matrix.cols());
        assert(static_cast<size_t>(column_coordinate_begin.row) < debug_score_matrix.rows());

        // We have two algorithms:
        // First copy the alignment
        std::ranges::copy(alignment_column | std::views::transform([] (auto && cell) { return cell.best_score(); }),
                          debug_score_matrix.begin() + column_coordinate_begin);

        if constexpr (traits_type::compute_sequence_alignment)
        {
            assert(static_cast<size_t>(column_coordinate_begin.col) < debug_trace_matrix.cols());
            assert(static_cast<size_t>(column_coordinate_begin.row) < debug_trace_matrix.rows());

            std::ranges::copy(alignment_column | std::views::transform([] (auto && cell) { return cell.best_trace(); }),
                              debug_trace_matrix.begin() + column_coordinate_begin);
        }
    }
};
} // namespace seqan3::detail
