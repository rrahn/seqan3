// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_storage_pair.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{
/*!\brief An alignment matrix that combines a score matrix with a trace matrix into a common interface.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_matrix_t The type of the score matrix; must model std::ranges::input_range (for more requirements see
 *                        details section).
 * \tparam trace_matrix_t The type of the trace matrix; must model std::ranges::input_range (for more requirements see
 *                        details section).
 *
 * \details
 *
 * To compute the actual, alignment an additional trace matrix is required. It stores the trace directions
 * to the previous cells for any computed cell within the alignment matrix. The following matrix combines a given score
 * with a trace matrix into a complete alignment matrix. The iterator increments over the columns of each underlying
 * matrix and combines them into a zip view, i.e. the combined alignment matrix is a range over ranges with the inner
 * range also modelling std::ranges::input_range.
 * To provide a uniform interface to the alignment algorithm, the returned zip view is further transformed such that
 * its reference type is a seqan3::detail::affine_cell_proxy.
 */
template <typename matrix_storage_first_t, typename matrix_storage_second_t>
class alignment_matrix_storage_pair
{
private:

    class iterator;
    class sentinel;

    //!\brief The underlying score matrix.
    matrix_storage_first_t first_storage{};
    //!\brief The underlying trace matrix.
    matrix_storage_second_t second_storage{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_matrix_storage_pair() = default; //!< Defaulted.
    alignment_matrix_storage_pair(alignment_matrix_storage_pair const &) = default; //!< Defaulted.
    alignment_matrix_storage_pair(alignment_matrix_storage_pair &&) = default; //!< Defaulted.
    alignment_matrix_storage_pair & operator=(alignment_matrix_storage_pair const &) = default; //!< Defaulted.
    alignment_matrix_storage_pair & operator=(alignment_matrix_storage_pair &&) = default; //!< Defaulted.
    ~alignment_matrix_storage_pair() = default; //!< Defaulted.
    //!\}

    /*!\brief Resizes the matrix.
     * \tparam column_index_t The column index type; must model std::integral.
     * \tparam row_index_t The row index type; must model std::integral.
     *
     * \param[in] column_count The number of columns for this matrix.
     * \param[in] row_count The number of rows for this matrix.
     * \param[in] initial_score The initial score used to initialise the score matrix.
     *
     * \details
     *
     * Resizes the underlying score and trace matrix to the given dimensions.
     *
     * ### Complexity
     *
     * The complexity depends on the resize complexity of the underlying score and trace matrix.
     *
     * ### Exception
     *
     * Strong exception guarantee. Might throw std::bad_alloc.
     */
    template <std::integral row_index_t, std::integral column_index_t, std::semiregular score_t>
    void resize(row_index_type<row_index_t> const row_count,
                column_index_type<column_index_t> const column_count,
                score_t const initial_score = score_t{})
    {
        score_matrix_t tmp_score_matrix{};
        tmp_score_matrix.resize(row_count, column_count, initial_score);

        trace_matrix_t tmp_trace_matrix{};
        tmp_trace_matrix.resize(row_count, column_count);

        score_matrix = std::move(tmp_score_matrix);
        trace_matrix = std::move(tmp_trace_matrix);
    }
    // When we zip it, how is the value type?
    // value_type = std::pair<range_value_t<score_col_t>, range_value_t<trace_col_t>>
    //  => std::pair<affine_cell<std::tuple<int, int, int, int>>, trace_directions>
    // reference = ranges::common_pair<range_reference_t<score_col_t>, range_reference_t<trace_col_t>>
    //  => ranges::common_pair<affine_cell<std::tuple<int &, int &, int &, int &>>, trace_directions &>
    // assignable from:
    // std::pair{affine_cell{0, 1, 2}, trace_directions::none};
    // vector over value type should be assignable from the range

    // But we do have a concrete type that is expected inside of the alignment.
    // This is the score + trace
    column_type column_at(size_t const column_id)
    {
        // We could make it a proxy type:
        // first = score_cell()
        // second = trace;
        return views::zip(score_matrix.column_at(column_id), trace_matrix.column_at(column_id));
    }

    /*!\brief Forwards the request to the contained trace matrix.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     *
     * \returns An iterator pointing to the underlying trace matrix.
     *
     * \throws std::out_of_range if the specified coordinate is out of range.
     */
    auto matrix_iterator_at(matrix_coordinate const & trace_begin) const
    {
        return trace_matrix.matrix_iterator_at(trace_begin);
    }
};

} // namespace seqan3::detail


// Design A) The trace matrix is never stand alone but always in combination with the score matrix.
// pro: No extra comparison for the trace matrix if we know already that they have the same size.
//      Extra data member needed for the score matrix anyway since the recursion is different.
//
// con:

// Design B) The trace matrix is a separate entity and is later combined.
// pro:      Combination is part of another abstraction and not of the trace matrix.
//
// con:      We need to specialise the score matrix type anyway.

