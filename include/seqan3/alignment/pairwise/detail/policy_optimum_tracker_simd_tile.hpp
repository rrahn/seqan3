// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_optimum_tracker.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <seqan3/std/ranges>

#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker_simd.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

template <simd_concept simd_index_t>
struct alignment_optimum_updater_greater_equal_global_alignment_simd_2
{
    //!\brief The column indices that need to match.
    simd_index_t column_index{};
    //!\brief The row indices that need to match.
    simd_index_t row_index{};

    typename simd_traits<simd_index_t>::mask_type track_this_tile{};

    /*!\brief The binary compare and update operation.
     * \tparam lhs_t The type of the left hand side; must model seqan3::tuple_like with a tuple size of 2.
     * \tparam rhs_t The type of the right hand side; must model seqan3::tuple_like with a tuple size of 2.
     *
     * \param[in,out] optimal_score_coordinate_pair The current optimum.
     * \param[in] current_cell_score_coordinate_pair The new score and coordinate to compare with.
     *
     * \details
     *
     * Requires that first value of the tuple represents the simd score and the second type the matrix coordinate over
     * the simd column and row index.
     */
    template <tuple_like lhs_t, tuple_like rhs_t>
    // //!\cond
    //     requires std::tuple_size_v<lhs_t> == 2 &&
    //              std::tuple_size_v<rhs_t> == 2 &&
    //              simd_concept<std::tuple_element_t<0, lhs_t>> &&
    //              simd_concept<std::tuple_element_t<0, rhs_t>>
    // //!\endcond
    void operator()(lhs_t && optimal_score_coordinate_pair,
                    rhs_t && current_cell_score_coordinate_pair) const
    {
        auto && [optimal_score, optimal_coordinate] = optimal_score_coordinate_pair;
        auto && [current_cell_score, current_cell_coordinate] = current_cell_score_coordinate_pair;

        auto mask = track_this_tile && (column_index == current_cell_coordinate.col) &&
                    (row_index == current_cell_coordinate.row);
        optimal_score = (mask) ? current_cell_score : optimal_score;
        optimal_coordinate.col = (mask) ? column_index : optimal_coordinate.col;
        optimal_coordinate.row = (mask) ? row_index : optimal_coordinate.row;
    }
};


/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_optimum_tracker
 */
template <typename alignment_configuration_t, std::semiregular binary_update_operation_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class policy_optimum_tracker_simd_tile :
    protected policy_optimum_tracker_simd<alignment_configuration_t, binary_update_operation_t>
{
protected:
    //!\brief The type of the base class.
    using base_policy_t = policy_optimum_tracker_simd<alignment_configuration_t, binary_update_operation_t>;

    // Import the configured score type.
    using typename base_policy_t::traits_type;
    using typename base_policy_t::score_type;
    using typename base_policy_t::matrix_coordinate_type;
    using typename base_policy_t::matrix_index_type;

    using original_score_type = typename traits_type::original_score_type;
    // We need to make this traits types.
    using cell_type = affine_cell_proxy<std::tuple<score_type, score_type, score_type>>;
    using cell_vector_type = std::vector<cell_type, aligned_allocator<cell_type, alignof(cell_type)>>;
    using tile_iterator_type = std::ranges::iterator_t<cell_vector_type>;


    using base_policy_t::optimal_score;
    using base_policy_t::optimal_coordinate;
    using base_policy_t::binary_update_operation;
    using base_policy_t::padding_offsets;

    tile_iterator_type tile_column_it{};
    tile_iterator_type tile_row_it{};

    // Mark the column index of a tile that contains an optimum to track.
    matrix_index_type tile_column_index{};
    // Mark the row index of a tile that contains an optimum to track.
    matrix_index_type tile_row_index{};

    std::array<original_score_type, traits_type::alignments_per_vector> score_column_offsets{};
    std::array<original_score_type, traits_type::alignments_per_vector> score_row_offsets{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_optimum_tracker_simd_tile() = default; //!< Defaulted.
    policy_optimum_tracker_simd_tile(policy_optimum_tracker_simd_tile const &) = default; //!< Defaulted.
    policy_optimum_tracker_simd_tile(policy_optimum_tracker_simd_tile &&) = default; //!< Defaulted.
    policy_optimum_tracker_simd_tile & operator=(policy_optimum_tracker_simd_tile const &) = default; //!< Defaulted.
    policy_optimum_tracker_simd_tile & operator=(policy_optimum_tracker_simd_tile &&) = default; //!< Defaulted.
    ~policy_optimum_tracker_simd_tile() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration (not used in this context).
     *
     * \details
     *
     * Resets the optimum on construction. Sets track_last_row_cell and track_last_column_cell to true in order
     * to check if this is a cell with a potential optimum.
     */
    policy_optimum_tracker_simd_tile(alignment_configuration_t const & config) : base_policy_t{config}
    {}
    //!\}

    template <typename cell_t>
    decltype(auto) track_last_row_cell(cell_t && cell, matrix_coordinate_type SEQAN3_DOXYGEN_ONLY(coordinate)) noexcept
    {
        // We store the cell before we call update.
        *tile_row_it = cell;
        ++tile_row_it;
        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_last_column_cell(cell_t && cell, matrix_coordinate_type SEQAN3_DOXYGEN_ONLY(coordinate)) noexcept
    {
        // We store the cell before we call update.
        *tile_column_it = cell;
        ++tile_column_it;
        return std::forward<cell_t>(cell);
    }

    template <typename cell_vector_t, typename tile_column_offsets_t>
    void track_tile_last_row(cell_vector_t & last_row,
                             matrix_coordinate_type const tile_coordinate,
                             original_score_type const row_offset,
                             tile_column_offsets_t const & tile_column_offsets)
    {
        // What happens now:
        // We get the coordinate.
        auto mask = tile_coordinate.col == tile_column_index && tile_coordinate.row == tile_row_index;
        // We only track the cells if this is true.
        binary_update_operation.track_this_tile = mask;

        // Now we go over each cell in the last row and track the optimum if this tile needs to be tracked and
        // if the column index corresponds to the local column index.
        // We set the row index equal to the stored row_index such that this is always true.

        for (auto && [cell, column_index] : views::zip(last_row, views::iota_simd<matrix_index_type>(0u, last_row.size())))
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), matrix_coordinate_type{row_index_type{binary_update_operation.row_index}, column_index_type{column_index}}));

        // Go over the mask and track the row offsets as well as the column offsets for every tracked element.
        for (size_t idx = 0; idx < simd_traits<matrix_index_type>::length; ++idx)
        {
            if (!mask[idx])
                continue;

            // Depends on the initialisation?
            score_row_offsets[idx] = row_offset;

            for (auto const & column_offset : tile_column_offsets)
                score_column_offsets[idx] += column_offset[idx];
        }
    }

    template <typename cell_vector_t, typename tile_column_offsets_t>
    void track_tile_last_column(cell_vector_t & last_column,
                                matrix_coordinate_type const tile_coordinate,
                                original_score_type const row_offset,
                                tile_column_offsets_t const & tile_column_offsets)
    {
        // What happens now:
        // We get the coordinate.
        auto mask = tile_coordinate.col == tile_column_index && tile_coordinate.row == tile_row_index;
        // We only track the cells if this is true.
        binary_update_operation.track_this_tile = mask;

        // Now we go over each cell in the last column and track the optimum if this tile needs to be tracked and
        // if the column index corresponds to the local column index.
        // We set the row index equal to the stored row_index such that this is always true.
        for (auto && [cell, row_index] : views::zip(last_column, views::iota_simd<matrix_index_type>(0u, last_column.size())))
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), matrix_coordinate_type{row_index_type{row_index}, column_index_type{binary_update_operation.column_index}}));

        // Go over the mask and track the row offsets as well as the column offsets for every tracked element.
        for (size_t idx = 0; idx < simd_traits<matrix_index_type>::length; ++idx)
        {
            if (!mask[idx])
                continue;

            score_row_offsets[idx] = row_offset;

            for (size_t i = 0; i < tile_coordinate.row[idx]; ++i)
                score_column_offsets[idx] += tile_column_offsets[i][idx];
        }
    }

    void set_tile_to_track(cell_vector_type & column,
                           cell_vector_type & row,
                           matrix_coordinate_type const simd_coordinate)
    {
        binary_update_operation.track_this_tile = (simd_coordinate.col == tile_column_index) &&
                                                  (simd_coordinate.row == tile_row_index);
        tile_column_it = column.begin();
        tile_row_it = row.begin();
    }

    auto row_offset(size_t const index)
    {
        return score_row_offsets[index]; // Who stores this?
    }

    auto column_offset(size_t const index)
    {
        return score_column_offsets[index]; // Who stores this?
    }

    template <std::ranges::input_range sequence1_collection_t, std::ranges::input_range sequence2_collection_t>
    void initialise_tracker(sequence1_collection_t & sequence1_collection,
                            sequence2_collection_t & sequence2_collection,
                            size_t const tile_size)
    {
        using scalar_index_t = uint16_t; // We set the maximum sequence size to uint16_t.
        using simd_uint16_t = simd_type_t<scalar_index_t>;  // The corresponding simd type.
        //
        scalar_index_t largest_sequence1_size{};
        scalar_index_t largest_sequence2_size{};
        alignas(alignof(simd_uint16_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence1_sizes{};
        alignas(alignof(simd_uint16_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence2_sizes{};

        // First, get all dimensions from the sequences and keep track of the maximal sizes.
        size_t index{};
        for (auto && [sequence1, sequence2] : views::zip(sequence1_collection, sequence2_collection))
        {
            sequence1_sizes[index] = static_cast<scalar_index_t>(std::ranges::distance(sequence1));
            sequence2_sizes[index] = static_cast<scalar_index_t>(std::ranges::distance(sequence2));
            largest_sequence1_size = std::max(largest_sequence1_size, sequence1_sizes[index]);
            largest_sequence2_size = std::max(largest_sequence2_size, sequence2_sizes[index]);
            ++index;
        }

        assert(index > 0); // We do not expect empty sequence collections here.

        // Second, get the offset for each individual end coordinate to project the cell to the last row or
        // column of the global alignment matrix. Choose the smallest distance, since this gives the correct offset
        // to the projected end cell.
        do
        {
            --index;

            assert(sequence1_sizes[index] <= largest_sequence1_size);
            assert(sequence2_sizes[index] <= largest_sequence2_size);

            padding_offsets[index] = std::min((largest_sequence1_size - sequence1_sizes[index]),
                                              (largest_sequence2_size - sequence2_sizes[index]));
            sequence1_sizes[index] += padding_offsets[index];
            sequence2_sizes[index] += padding_offsets[index];
        } while (index > 0);

        // Load the target coordinate indices from the respective arrays.
        // Either we do computation sequential or we compute them in 2 vector types.

        simd_uint16_t column_index_lo = simd::load<simd_uint16_t>(sequence1_sizes.data());
        simd_uint16_t column_index_hi = simd::load<simd_uint16_t>(sequence1_sizes.data() + simd_traits<simd_uint16_t>::length);

        simd_uint16_t row_index_lo = simd::load<simd_uint16_t>(sequence2_sizes.data());
        simd_uint16_t row_index_hi = simd::load<simd_uint16_t>(sequence2_sizes.data() + simd_traits<simd_uint16_t>::length);

        simd_uint16_t simd_tile_size = simd::fill<simd_uint16_t>(tile_size);
        simd_uint16_t simd_base = simd_tile_size - simd::fill<simd_uint16_t>(1);
        simd_uint16_t mask_one = simd::fill<simd_uint16_t>(1);

        // Now transform the computed column and row indices to the relative tile indices:
        // First compute the relative tile index. We assume that the relative tile index fits into an uint8_t type.
        // Second pack the 16 bit of hi and lo index into a simd vector with 8 bit packed integers. The pack_epi16
        // instruction puts the first 8x16 bit integers of `a` into the first 128 bit followed by the first 8x16 bit
        // integers of `b` and then again with the remaining 8x16 bit integers of `a` followed by the ones from `b`.
        // Accordingly, the bits in the range from [127...64] must be swapped with the bits in the range [191...128].
        // Third shuffle the remaining vector, i.e: [63...0] -> 0, [127...64] -> 2, [191...128] -> 1, [255...192] -> 3
        // The shuffle mask in binary format: 0b11'01'10'00.

        simd_uint16_t tile_column_index_hi = (column_index_hi + simd_base) % simd_tile_size + mask_one;
        simd_uint16_t tile_column_index_lo = (column_index_lo + simd_base) % simd_tile_size + mask_one;

        // debug_stream << "absolute: " << column_index_hi << " relative: " << tile_column_index_hi;

        binary_update_operation.column_index = reinterpret_cast<matrix_index_type>(
            _mm256_permute4x64_epi64(_mm256_packus_epi16(reinterpret_cast<__m256i &>(tile_column_index_lo),
                                                         reinterpret_cast<__m256i &>(tile_column_index_hi)),
                                     0b11011000));


        simd_uint16_t tile_row_index_hi = (row_index_hi + simd_base) % simd_tile_size + mask_one;
        simd_uint16_t tile_row_index_lo = (row_index_lo + simd_base) % simd_tile_size + mask_one;

        // debug_stream << "absolute: " << row_index_hi << " relative: " << tile_row_index_hi;

        binary_update_operation.row_index = reinterpret_cast<matrix_index_type>(
            _mm256_permute4x64_epi64(_mm256_packus_epi16(reinterpret_cast<__m256i &>(tile_row_index_lo),
                                                         reinterpret_cast<__m256i &>(tile_row_index_hi)),
                                     0b11011000));

        // For now we assume that 256*tile_size is enough to index the tiles.
        // The index is also stored as a 8 bit packed simd vector.
        // We proceed like above.
        // First calculate the tile index and then store the results in the respective vector.

        tile_column_index_hi = (column_index_hi + simd_base) / simd_tile_size - mask_one;
        tile_column_index_lo = (column_index_lo + simd_base) / simd_tile_size - mask_one;

        // debug_stream << "absolute: " << column_index_hi << " relative: " << tile_column_index_hi;

        tile_column_index = reinterpret_cast<matrix_index_type>(
            _mm256_permute4x64_epi64(_mm256_packus_epi16(reinterpret_cast<__m256i &>(tile_column_index_lo),
                                                         reinterpret_cast<__m256i &>(tile_column_index_hi)),
                                     0b11011000));

        tile_row_index_hi = (row_index_hi + simd_base) / simd_tile_size - mask_one;
        tile_row_index_lo = (row_index_lo + simd_base) / simd_tile_size - mask_one;

        // debug_stream << "absolute: " << row_index_hi << " relative: " << tile_row_index_hi;

        tile_row_index = reinterpret_cast<matrix_index_type>(
            _mm256_permute4x64_epi64(_mm256_packus_epi16(reinterpret_cast<__m256i &>(tile_row_index_lo),
                                                         reinterpret_cast<__m256i &>(tile_row_index_hi)),
                                     0b11011000));

        // Initialise the offset fields:
        score_column_offsets.fill(original_score_type{});
        score_row_offsets.fill(original_score_type{});
    }

};
} // namespace seqan3::detail
