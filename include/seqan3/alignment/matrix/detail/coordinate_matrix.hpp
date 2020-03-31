// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::coordinate_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>

namespace seqan3::detail
{

/*!\brief A function object that converts a column index and a row index to a seqan3::detail::matrix_coordinate.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This helper function object is used for the seqan3::detail::coordinate_matrix_iterator to convert the a column index
 * and row index into a seqan3::detail::matrix_coordinate. The seqan3::detail::coordinate_matrix_iterator iterates
 * over the columns of the underlying coordinate matrix and returns a transformed iota range over the row indices of
 * this column. The column index is stored as member inside of this function object since every column index is
 * associated with many row indices.
 */
struct convert_to_matrix_coordinate
{
    //!\brief The index of the represented column.
    size_t column_index{0};

    /*!\brief The conversion operator
     *
     * \param[in] row_index The index of the represented row.
     *
     * \returns seqan3::detail::matrix_coordinate for the current column and row index.
     */
    auto operator()(size_t const row_index) noexcept
    {
        return matrix_coordinate{row_index_type{row_index}, column_index_type{column_index}};
    }
};

/*!\brief The iterator for the seqan3::detail::coordinate_matrix.
 * \ingroup alignment_matrix
 * \implements std::forward_iterator
 *
 * \details
 *
 * Iterates over the columns of the underlying seqan3::detail::coordinate_matrix. The iterator returns a transformed
 * std::ranges::views::iota view over the row indices. In the transformation every row index is converted to a
 * seqan3::detail::matrix_coordinate using the seqan3::detail::convert_to_matrix_coordinate function object.
 */
class coordinate_matrix_iterator
{
private:

    //!\brief The currently represented column index.
    size_t m_column_index{0};
    //!\brief The row index marking the end of the rows (size of one column).
    size_t m_end_row_index{0};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::views::iota(static_cast<size_t>(0), static_cast<size_t>(1))
                              | std::views::transform(convert_to_matrix_coordinate{m_column_index}));
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::forward_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    coordinate_matrix_iterator() noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator(coordinate_matrix_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator(coordinate_matrix_iterator &&) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator & operator=(coordinate_matrix_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator & operator=(coordinate_matrix_iterator &&) noexcept = default; //!< Defaulted.
    ~coordinate_matrix_iterator() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the current column index and the row index marking the end
     *        of the rows (size of one column).
     *
     * \param[in] column_index \copybrief seqan3::detail::coordinate_matrix_iterator::m_column_index
     * \param[in] end_row_index \copybrief seqan3::detail::coordinate_matrix_iterator::m_end_row_index
     */
    explicit coordinate_matrix_iterator(column_index_type<size_t> column_index,
                                        row_index_type<size_t> end_row_index) noexcept :
        m_column_index{column_index.get()},
        m_end_row_index{end_row_index.get()}
    {}
    //!\}

    /*!\name Access operators
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    auto operator*() const
    {
        return std::views::iota(static_cast<size_t>(0), m_end_row_index)
             | std::views::transform(convert_to_matrix_coordinate{m_column_index});
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    coordinate_matrix_iterator & operator++()
    {
        ++m_column_index;
        return *this;
    }

    //!\brief Increments the iterator to the next column and returns the iterator pointing to the previous column.
    coordinate_matrix_iterator operator++(int)
    {
        coordinate_matrix_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(coordinate_matrix_iterator const & lhs, coordinate_matrix_iterator const & rhs)
    {
        return lhs.m_column_index == rhs.m_column_index;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(coordinate_matrix_iterator const & lhs, coordinate_matrix_iterator const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}
};

/*!\brief A matrix over coordinates.
 * \ingroup alignment_matrix
 * \implements std::ranges::forward_range
 *
 * \details
 *
 * In the alignment algorithm an additional matrix is needed to represent the coordinates of the processed cells.
 * This matrix is cheap as it stores only the dimensions of the matrix and does not allocate any memory for the
 * coordinates. It uses the seqan3::detail::matrix_coordinate_iterator to iterate
 * over the columns of this virtual matrix. When the seqan3::detail::matrix_coordinate_iterator is dereferenced it
 * returns an on-the-fly constructed range ([prvalue](https://en.cppreference.com/w/cpp/language/value_category))
 * representing the current coordinate column. The value type of this range is seqan3::detail::matrix_coordinate.
 */
class coordinate_matrix
{
private:

    //!\brief The column index marking the end of columns (size of one row).
    column_index_type<size_t> m_end_column_index{};
    //!\brief The row index marking the end of rows (size of one column).
    row_index_type<size_t> m_end_row_index{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    coordinate_matrix() = default; //!< Defaulted.
    coordinate_matrix(coordinate_matrix const &) = default; //!< Defaulted.
    coordinate_matrix(coordinate_matrix &&) = default; //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix const &) = default; //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix &&) = default; //!< Defaulted.
    ~coordinate_matrix() = default; //!< Defaulted.

    /*!\brief Resets the coordinate matrix with the given end column index and end row index representing the dimensions
     *        of the reset matrix.
     *
     * \param[in] end_column_index \copybrief seqan3::detail::coordinate_matrix::m_end_column_index
     * \param[in] end_row_index \copybrief seqan3::detail::coordinate_matrix::m_end_row_index
     */
    void reset_matrix(column_index_type<size_t> const end_column_index, row_index_type<size_t> const end_row_index)
    {
        m_end_column_index = end_column_index;
        m_end_row_index = end_row_index;
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column of the matrix.
    coordinate_matrix_iterator begin() const noexcept
    {
        return coordinate_matrix_iterator{column_index_type<size_t>{0u}, m_end_row_index};
    }

    //!\brief Returns the iterator pointing to the end column of the matrix.
    coordinate_matrix_iterator end() const noexcept
    {
        return coordinate_matrix_iterator{m_end_column_index, m_end_row_index};
    }
    //!\}
};

} // namespace seqan3::detail
