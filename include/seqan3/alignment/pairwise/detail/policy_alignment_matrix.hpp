// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <tuple>
#include <seqan3/std/concepts>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/detail/combined_score_and_trace_matrix.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_simd_adapter_iterator.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief A policy that provides a common interface to acquire the correct alignment matrices.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_matrix_t The type of the alignment matrix [see requirements below].
 * \tparam coordinate_matrix_t The type of the coordinate matrix [see requirements below].
 *
 * \details
 *
 * The alignment matrix must be a matrix type that is compatible with the configured alignment algorithm. It must offer
 * a resize member function that takes a seqan3::detail::column_index_type and seqan3::detail::row_index_type and an
 * additional parameter to initialise the allocated matrix memory.
 */
template <typename alignment_matrix_t, typename coordinate_matrix_t>
//!\cond
    requires requires (alignment_matrix_t & alignment_matrix,
                       coordinate_matrix_t & index_matrix,
                       typename alignment_matrix_t::score_type const initial_score)
             {
                 { alignment_matrix.resize(column_index_type{size_t{}}, row_index_type{size_t{}}, initial_score) };
                 { index_matrix.resize(column_index_type{size_t{}}, row_index_type{size_t{}}) };
             }
//!\endcond
class policy_alignment_matrix
{
protected:
    // ----------------------------------------------------------------------------
    // type definition
    // ----------------------------------------------------------------------------

    //!\brief The underlying score type of the alignment matrix.
    using score_type = typename alignment_matrix_t::score_type;

    // ----------------------------------------------------------------------------
    // static member
    // ----------------------------------------------------------------------------

    //!\brief Flag indicating whether the trace was computed.
    static constexpr bool with_trace = is_type_specialisation_of_v<alignment_matrix_t, combined_score_and_trace_matrix>;
    //!\brief Flag indicating whether the alignment is computed in vectorised mode.
    static constexpr bool is_vectorised = simd_concept<score_type>;

    // ----------------------------------------------------------------------------
    // non-static member
    // ----------------------------------------------------------------------------

    //!\brief The selected lower diagonal.
    int32_t lower_diagonal{};
    //!\brief The selected upper diagonal.
    int32_t upper_diagonal{};
    //!\brief The largest column index where the alignment matrix (banded or not) still intersects in the first row.
    int32_t largest_valid_column_index{};
    //!\brief A flag indicating whether the final gaps in the last column are free.
    bool last_column_is_free{};
    //!\brief A flag indicating whether the final gaps in the last row are free.
    bool last_row_is_free{};
    //!\brief A flag indicating global alignment.
    bool is_global{};
    //!\brief A flag indicating banded alignment.
    bool is_banded{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_matrix() = default; //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix const &) = default; //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix &&) = default; //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix const &) = default; //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix &&) = default; //!< Defaulted.
    ~policy_alignment_matrix() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \tparam alignment_configuration_t The type of the alignment configuration; must be an instance of
     *                                   seqan3::configuration.
     *
     * \param[in] config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the members for the lower and upper diagonal. These members are only used if the banded alignment
     * is computed.
     *
     * \throws seqan3::invalid_alignment_configuration if the given band settings are invalid.
     */
    template <typename alignment_configuration_t>
    //!\cond
        requires (is_type_specialisation_of_v<alignment_configuration_t, configuration>)
    //!\endcond
    policy_alignment_matrix(alignment_configuration_t const & config)
    {
        using seqan3::get;

        is_global = config.template exists<align_cfg::method_global>();
        is_banded = config.template exists<align_cfg::band_fixed_size>();

        auto band = config.get_or(seqan3::align_cfg::band_fixed_size{});

        lower_diagonal = band.lower_diagonal;
        upper_diagonal = band.upper_diagonal;

        bool invalid_band = upper_diagonal < lower_diagonal;
        std::string error_cause = (invalid_band) ? " The upper diagonal is smaller than the lower diagonal." : "";

        if (is_global)
        {
            auto method_global_config = config.get_or(align_cfg::method_global{});

            bool first_row_is_free = method_global_config.free_end_gaps_sequence1_leading;
            bool first_column_is_free = method_global_config.free_end_gaps_sequence2_leading;

            last_row_is_free = method_global_config.free_end_gaps_sequence1_trailing;
            last_column_is_free = method_global_config.free_end_gaps_sequence2_trailing;
            // band starts in first column without free gaps or band starts in first row without free gaps.
            invalid_band |= (upper_diagonal < 0 && !first_column_is_free) || (lower_diagonal > 0 && !first_row_is_free);
            error_cause += " The band starts in a region without free gaps.";
        }

        if (invalid_band)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":" +
                                                  std::to_string(upper_diagonal) + "] cannot be used with the current "
                                                  "alignment configuration:" + error_cause};
    }
    //!\}

    /*!\brief Acquires a new thread local alignment and index matrix for the given sequence sizes.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     * \param[in] initial_score The initial score used for the acquired alignment matrix.
     *
     * \returns A std::tuple storing lvalue references to the thread local alignment and index matrix.
     *
     * \details
     *
     * Acquires a thread local alignment and index matrix. Initialises the matrices with the given
     * sequence sizes and the initial score value. In the banded alignment, the alignment matrix is reduced to
     * the column count times the band size.
     *
     * ### Exception
     *
     * Might throw std::bad_alloc if the requested matrix size exceeds the available memory or
     * seqan3::invalid_alignment_configuration if the band does not allow a valid computation of the
     * configured alignment.
     *
     * \throws std::bad_alloc or seqan3::invalid_alignment_configuration
     */
    auto acquire_matrices(size_t const sequence1_size,
                          size_t const sequence2_size,
                          score_type initial_score = score_type{})
    {
        assert(sequence1_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        assert(sequence2_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));

        largest_valid_column_index = std::clamp<int32_t>(upper_diagonal, 0, sequence1_size);

        if (is_banded)
            check_valid_band_configuration(sequence1_size, sequence2_size);

        static thread_local alignment_matrix_t alignment_matrix{};
        static thread_local coordinate_matrix_t coordinate_matrix{};

        // Increase dimension by one for the initialisation of the matrix.
        size_t const column_count = sequence1_size + 1;
        size_t row_count = sequence2_size + 1;

        coordinate_matrix.resize(column_index_type{column_count}, row_index_type{row_count});

        if (is_banded)
        {
            assert(upper_diagonal - lower_diagonal + 1 > 0); // Band size is a positive integer.
            // Allocate one more cell to compute the last cell of the band with standard recursion function.
            row_count = std::min<int64_t>(upper_diagonal - lower_diagonal + 2, row_count);
        }

        alignment_matrix.resize(column_index_type{column_count}, row_index_type{row_count}, initial_score);

        return std::tie(alignment_matrix, coordinate_matrix);
    }

    /*!\brief Checks whether the band is valid for the given sequence sizes.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     *
     * \throws seqan3::invalid_alignment_configuration if the band is invalid for the given sequence sizes and the
     *         alignment configuration.
     */
    void check_valid_band_configuration(size_t const sequence1_size, size_t const sequence2_size) const
    {
        bool const upper_diagonal_ends_before_last_cell = (upper_diagonal + sequence2_size) < sequence1_size;
        bool const lower_diagonal_ends_behind_last_cell = (-lower_diagonal + sequence1_size) < sequence2_size;

        bool invalid_band = false;
        std::string error_cause{};

        if (is_global)
        {
            // band ends in last column without free gaps or band ends in last row without free gaps.
            invalid_band |= (lower_diagonal_ends_behind_last_cell && !last_column_is_free) ||
                            (upper_diagonal_ends_before_last_cell && !last_row_is_free);
            error_cause = "The band ends in a region without free gaps.";
        }

        if (invalid_band)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":" +
                                                  std::to_string(upper_diagonal) + "] cannot be used with the current "
                                                  "alignment configuration: " + error_cause};
    }

    /*!\brief Returns the trace path starting at the given coordinate and for the respective simd position
     *
     * \param[in] alignment_matrix The alignment matrix to get the trace path from.
     * \param[in] coordinate The alignment matrix coordinate where the trace path starts.
     * \param[in] simd_position The simd vector position associated with the respective alignment matrix in vectorised
     *                          alignment.
     *
     * \details
     *
     * Creates a lazy range over the trace path starting at the given matrix coordinate. The end of the trace is
     * defined by the first coordinate who's trace value is equal to seqan3::detail::trace_directions::none.
     * In the banded alignment the row index of the given coordinate is adapted to map to the correct position within
     * the trace matrix, since only a fraction of the original memory was allocated.
     *
     * ### Exception
     *
     * Might throw seqan3::invalid_alignment_configuration if this function is invoked even if the computation of the
     * trace was not configured by the user.
     * Otherwise, might throw implementation defined exceptions.
     * If an exception is thrown the state of the policy is not changed (strong exception guarantee).
     */
    auto trace_path_starting_at(alignment_matrix_t const & alignment_matrix,
                                matrix_coordinate coordinate,
                                [[maybe_unused]] size_t const simd_position)
    {
        if constexpr (with_trace)
        {
            if (is_banded)  // update the coordinate's row index if the alignment was computed with a band.
            {
                int32_t column_index = coordinate.col;
                coordinate.row -= (column_index > largest_valid_column_index) *
                                  (column_index - largest_valid_column_index);
            }

            // Helper function to make the actual trace path for the given trace matrix iterator.
            auto make_trace_path = [&] (auto trace_matrix_it)
            {
                using trace_matrix_iterator_t = decltype(trace_matrix_it);
                using trace_path_iterator_t = trace_iterator<trace_matrix_iterator_t>;
                using trace_path_t = std::ranges::subrange<trace_path_iterator_t, std::default_sentinel_t>;

                return trace_path_t{trace_path_iterator_t{trace_matrix_it,
                                                          column_index_type{largest_valid_column_index}},
                                    std::default_sentinel};
            };

            auto trace_matrix_it = alignment_matrix.matrix_iterator_at(coordinate);

            if constexpr (is_vectorised)
                return make_trace_path(trace_matrix_simd_adapter_iterator{trace_matrix_it, simd_position});
            else
                return make_trace_path(trace_matrix_it);
        }
        else
        {
            static_assert(with_trace,
                          "Trying to get the trace path from an alignment without computed trace information.");
            return std::array<trace_directions, 0>{};
        }
    }
};
} // namespace seqan3::detail
