// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_matrix_full.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <seqan3/std/span>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

/*!\brief Trace matrix for the pairwise alignment using the full trace matrix.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam trace_t The type of the trace; must be the same as seqan3::detail::trace_directions.
 *
 * \details
 *
 * In the default trace back implementation we allocate the entire matrix using one byte per cell to store the
 * seqan3::detail::trace_directions.
 *
 * ### Range interface
 *
 * The matrix offers an input range interface over the columns of the matrix. Dereferencing the iterator will return
 * another range which represents the actual trace column in memory. The returned range is a
 * seqan3::views::zip view over the current column referencing the best trace, as well as the horizontal and vertical
 * trace column.
 */
template <typename trace_t>
class trace_matrix_full
{
private:
    //!\brief The type to store the complete trace matrix.
    using matrix_t = two_dimensional_matrix<trace_t,
                                            aligned_allocator<trace_t, sizeof(trace_t)>,
                                            matrix_major_order::column>;


    //!\brief The full trace matrix.
    matrix_t complete_matrix{};
    //!\brief The row count.
    size_t row_count{};

public:

    using column_type = std::span<trace_t>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    trace_matrix_full() = default; //!< Defaulted.
    trace_matrix_full(trace_matrix_full const &) = default; //!< Defaulted.
    trace_matrix_full(trace_matrix_full &&) = default; //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full const &) = default; //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full &&) = default; //!< Defaulted.
    ~trace_matrix_full() = default; //!< Defaulted.
    //!\}

    /*!\brief Resizes the matrix.
     * \tparam column_index_t The column index type; must model std::integral.
     * \tparam row_index_t The row index type; must model std::integral.
     *
     * \param[in] column_count The number of columns for this matrix.
     * \param[in] row_count The number of rows for this matrix.
     *
     * \details
     *
     * Resizes the entire trace matrix storing the best trace path and the horizontal trace column.
     * Note the trace matrix requires the number of columns and rows to be one bigger than the size of sequence1,
     * respectively sequence2 for the initialisation of the matrix.
     * Reallocation happens only if the new column size exceeds the current capacity of the underlying trace matrix.
     *
     * ### Complexity
     *
     * In worst case `column_count` times `row_count` memory is allocated.
     *
     * ### Exception
     *
     * Basic exception guarantee. Might throw std::bad_alloc on resizing the internal matrices.
     */
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(row_index_type<row_index_t> const row_count,
                column_index_type<column_index_t> const column_count)
    {
        this->row_count = row_count.get();
        complete_matrix.resize(number_rows{this->row_count}, number_cols{column_count.get()});
    }

    /*!\brief Returns an iterator pointing to the underlying trace matrix at the given matrix coordinate.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     * \returns An iterator pointing to the underlying trace matrix.
     * \throws std::out_of_range if the specified coordinate is out of range.
     */
    auto matrix_iterator_at(matrix_coordinate const & trace_begin) const
    {
        using matrix_iter_t = std::ranges::iterator_t<matrix_t const>;

        if (trace_begin.row >= row_count || trace_begin.col >= column_count)
            throw std::out_of_range{"The given coordinate exceeds the matrix in vertical or horizontal direction."};

        return matrix_iter_t{complete_matrix.begin() + matrix_offset{trace_begin}};
    }

    column_type column_at(size_t const column_id)
    {
        auto it = complete_matrix.begin() + matrix_offset{matrix_index{row_index_type{0}, column_index_type{column_id}}};
        return std::span<trace_t>{it, row_count};
    }
};

} // namespace seqan3::detail
