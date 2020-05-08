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
#include <seqan3/std/type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix_simd.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

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
    size_t column_index{};
    //!\brief The row index marking the end of the rows (size of one column).
    size_t end_row_index{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::views::iota(static_cast<size_t>(0), static_cast<size_t>(1))
                              | std::views::transform(convert_to_matrix_coordinate<size_t>{column_index}));
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
     * \param[in] column_index \copybrief seqan3::detail::coordinate_matrix_iterator::column_index
     * \param[in] end_row_index \copybrief seqan3::detail::coordinate_matrix_iterator::end_row_index
     */
    explicit coordinate_matrix_iterator(size_t column_index, size_t end_row_index) noexcept :
        column_index{column_index},
        end_row_index{end_row_index}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    auto operator*() const
    {
        return std::views::iota(static_cast<size_t>(0), end_row_index)
             | std::views::transform(convert_to_matrix_coordinate<size_t>{column_index});
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    coordinate_matrix_iterator & operator++()
    {
        ++column_index;
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
        return lhs.column_index == rhs.column_index;
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
 * \tparam index_t The underlying matrix index type; must mode std::integral or seqan3::simd::simd_concept.
 *
 * \details
 *
 * In the alignment algorithm an additional matrix is needed to represent the coordinates of the processed cells.
 * This matrix is cheap as it stores only the dimensions of the matrix and does not allocate any memory for the
 * coordinates. It uses the seqan3::detail::matrix_coordinate_iterator to iterate
 * over the columns of this virtual matrix. When the seqan3::detail::matrix_coordinate_iterator is dereferenced it
 * returns an on-the-fly constructed range ([prvalue](https://en.cppreference.com/w/cpp/language/value_category))
 * representing the current coordinate column. The reference type of this range is seqan3::detail::matrix_coordinate.
 *
 * ### Simd mode
 *
 * If the `index_t` is a simd vector the matrix implements a vectorised coordinate matrix instead. In this case the
 * matrix uses a seqan3::detail::matrix_coordinate_simd_iterator to iterate over the vectorised matrix. When the
 * seqan3::detail::matrix_coordinate_iterator is dereferenced it returns an on-the-fly constructed range
 * ([prvalue](https://en.cppreference.com/w/cpp/language/value_category)) representing the current simd coordinate
 * column. The reference type of this range is seqan3::detail::simd_matrix_coordinate.
 */
template <typename index_t>
//!\cond
    requires std::integral<index_t> || simd_concept<index_t>
//!\endcond
class coordinate_matrix
{
private:
    //!\brief Helper type alias to conditionally extract the scalar type from the seqan3::simd::simd_traits type.
    template <typename lazy_simd_traits_t>
    using lazy_scalar_type = typename instantiate_t<lazy_simd_traits_t>::scalar_type;

    //!\brief The internal size type which depends on `index_t` being a simd vector or a scalar type.
    using size_type = lazy_conditional_t<simd_concept<index_t>,
                                         lazy<lazy_scalar_type, lazy<simd_traits, index_t>>,
                                         size_t>;

    //!\brief The iterator type which depends on `index_t` being a simd vector or a scalar type.
    using iterator_type = lazy_conditional_t<simd_concept<index_t>,
                                             lazy<coordinate_matrix_simd_iterator, index_t>,
                                             coordinate_matrix_iterator>;

    //!\brief The column index marking the end of columns (size of one row).
    size_type end_column_index{};
    //!\brief The row index marking the end of rows (size of one column).
    size_type end_row_index{};

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

    /*!\brief Resets the coordinate matrix with the given end column index and end row index representing the new
     *        dimensions of the matrix.
     *
     * \param[in] end_column_index \copybrief seqan3::detail::coordinate_matrix::end_column_index
     * \param[in] end_row_index \copybrief seqan3::detail::coordinate_matrix::end_row_index
     */
    template <std::integral column_size_t, std::integral row_size_t>
    void reset_matrix(column_index_type<column_size_t> const end_column_index,
                      row_index_type<row_size_t> const end_row_index)
    {
        this->end_column_index = static_cast<size_type>(end_column_index.get());
        this->end_row_index = static_cast<size_type>(end_row_index.get());
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column of the matrix.
    iterator_type begin() const noexcept
    {
        return iterator_type{static_cast<size_type>(0), end_row_index};
    }

    //!\brief Returns the iterator pointing to the end column of the matrix.
    iterator_type end() const noexcept
    {
        return iterator_type{end_column_index, end_row_index};
    }
    //!\}
};

} // namespace seqan3::detail
