// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_recursion_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion_simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the vectorised alignment algorithm using affine gap costs.
 * \ingroup pairwise_alignment
 *
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename alignment_configuration_t>
class policy_affine_gap_recursion_simd_tile : protected policy_affine_gap_recursion_simd<alignment_configuration_t>
{
protected:
    //!\brief The type of the base class.
    using base_policy_t = policy_affine_gap_recursion_simd<alignment_configuration_t>;

    // Import the score type from the base.
    using typename base_policy_t::traits_type;
    using typename base_policy_t::score_type;
    using typename base_policy_t::affine_cell_type;
    //!\brief The original score type.
    using original_score_type = typename traits_type::original_score_type;
    //!\brief The matrix coordinate type.
    using matrix_coordinate_type = typename traits_type::matrix_coordinate_type;
    //!\brief The type of the row and column to reset.
    using cell_vector_type = std::vector<affine_cell_type, aligned_allocator<affine_cell_type, alignof(affine_cell_type)>>;
    //!\brief The iterator over the cell vector type.
    using cell_vector_iterator_type = std::ranges::iterator_t<cell_vector_type>;

    // Import parameters.
    using base_policy_t::gap_extension_score;
    using base_policy_t::gap_open_score;

    // Define them outside.
    // But where do we know about the types here?
    cell_vector_iterator_type tile_row_it{};
    cell_vector_iterator_type tile_column_it{};
    score_type tile_row_offset{};
    score_type tile_column_offset{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_recursion_simd_tile() = default; //!< Defaulted.
    policy_affine_gap_recursion_simd_tile(policy_affine_gap_recursion_simd_tile const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd_tile(policy_affine_gap_recursion_simd_tile &&) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd_tile & operator=(policy_affine_gap_recursion_simd_tile const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd_tile & operator=(policy_affine_gap_recursion_simd_tile &&) = default; //!< Defaulted.
    ~policy_affine_gap_recursion_simd_tile() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    explicit policy_affine_gap_recursion_simd_tile(alignment_configuration_t const & config) : base_policy_t{config}
    {}
    //!\}

    //!\brief overload
    affine_cell_type initialise_origin_cell() const noexcept
    {
        return apply_offset(*tile_column_it, tile_column_offset);
    }

    //!\brief overload
    template <typename affine_cell_t>
    //!\cond
        requires is_type_specialisation_of_v<affine_cell_t, affine_cell_proxy>
    //!\endcond
    affine_cell_type initialise_first_column_cell(affine_cell_t SEQAN3_DOXYGEN_ONLY(previous_cell)) noexcept
    {
        return apply_offset(*++tile_column_it, tile_column_offset);
    }

    //!\brief overload
    template <typename affine_cell_t>
    //!\cond
        requires is_type_specialisation_of_v<affine_cell_t, affine_cell_proxy>
    //!\endcond
    affine_cell_type initialise_first_row_cell(affine_cell_t SEQAN3_DOXYGEN_ONLY(previous_cell)) noexcept
    {
        return apply_offset(*++tile_row_it, tile_row_offset);
    }

    size_t determine_tile_size(int8_t const match_score, int8_t const mismatch_score) const
    {
        size_t candidate1 = std::numeric_limits<int8_t>::max() / (match_score - mismatch_score);
        size_t candidate2 = (std::numeric_limits<int8_t>::max() + gap_open_score[0]) /
                            (match_score - gap_extension_score[0]);
        return std::max(candidate1, candidate2);
    }

    void initialise_first_tile(cell_vector_type & column, cell_vector_type & row)
    {
        tile_row_offset = score_type{};
        tile_column_offset = score_type{};

        auto column_it = column.begin();
        *column_it = base_policy_t::initialise_origin_cell();
        ++column_it;
        for (; column_it != column.end(); ++column_it)
            *column_it = base_policy_t::initialise_first_column_cell(*std::ranges::prev(column_it));

        auto row_it = row.begin();
        *row_it = column.front();
        ++row_it;
        for (; row_it != row.end(); ++row_it)
            *row_it = base_policy_t::initialise_first_row_cell(*std::ranges::prev(row_it));

        set_tile(column, row);
    }

    void initialise_first_column_tile(cell_vector_type & column, cell_vector_type & row)
    {
        tile_row_offset = row.front().optimal_score();
        tile_column_offset = tile_row_offset;

        // Update the column.
        auto column_it = column.begin();
        *column_it = row.front();
        ++column_it;
        for (; column_it != column.end(); ++column_it)
            *column_it = base_policy_t::initialise_first_column_cell(*std::ranges::prev(column_it));

        set_tile(column, row);
    }

    void initialise_first_row_tile(cell_vector_type & column, cell_vector_type & row)
    {
        tile_column_offset = column.front().optimal_score();
        tile_row_offset = tile_column_offset;

        // Update the row.
        auto row_it = row.begin();
        *row_it = column.front();
        ++row_it;
        for (; row_it != row.end(); ++row_it)
            *row_it = base_policy_t::initialise_first_row_cell(*std::ranges::prev(row_it));

        set_tile(column, row);
    }

    void initialise_inner_tile(cell_vector_type & column, cell_vector_type & row)
    {
        tile_column_offset = column.front().optimal_score();
        tile_row_offset = row.front().optimal_score();

        set_tile(column, row);
    }

    // TODO: Move tile size into config.
    original_score_type compute_column_offset(matrix_coordinate_type const tile_coordinate, size_t const tile_size)
    {
        uint8_t current_column_index = tile_coordinate.col[0];
        // Compute the gap extension offset.
        original_score_type gap_extension_offset = current_column_index * tile_size *
                                                   original_score_type{gap_extension_score[0]};

        // Only add the gap open score if current column is > 0.
        return gap_extension_offset +
               (original_score_type{gap_open_score[0]} - original_score_type{gap_extension_score[0]}) *
               (current_column_index != 0);
    }

private:

    void set_tile(cell_vector_type & column, cell_vector_type & row)
    {
        tile_row_it = row.begin();
        tile_column_it = column.begin();
    }

    template <typename cell_t>
    affine_cell_type apply_offset(cell_t && cell, score_type const offset) const
    {
        return std::apply([&] (auto && ... cell_scores)
               {
                    return affine_cell_type{(cell_scores - offset)...};
               }, std::forward<cell_t>(cell));
    }
};
} // namespace seqan3::detail
