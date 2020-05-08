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

#include <seqan3/std/concepts>
#include <numeric>
#include <seqan3/std/ranges>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/view_iota_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

/*!\brief A function object that converts a column index and a row index to a seqan3::detail::matrix_coordinate.
 * \ingroup alignment_matrix
 *
 * \tparam index_t The type of the column and row index; must model std::integral or seqan3::simd::simd_concept.
 *
 * \details
 *
 * This helper function object is used for the seqan3::detail::coordinate_matrix_iterator and
 * seqan3::detail::coordinate_matrix_simd_iterator to convert the column and row index into a
 * seqan3::detail::matrix_index.
 */
template <typename index_t>
//!\cond
    requires std::integral<index_t> || simd_concept<index_t>
//!\endcond
struct convert_to_matrix_coordinate
{
    //!\brief The index of the represented column.
    index_t column_index{};

    /*!\brief The conversion operator
     *
     * \param[in] row_index The index of the represented row.
     *
     * \returns seqan3::detail::matrix_index for the current column and row index.
     */
    auto operator()(index_t const row_index) noexcept
    {
        return matrix_index{row_index_type{row_index}, column_index_type{column_index}};
    }
};

/*!\brief The iterator for the vectorised seqan3::detail::coordinate_matrix.
 * \ingroup alignment_matrix
 * \implements std::forward_iterator
 *
 * \tparam index_t The type of the index; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Iterates over the columns of the underlying seqan3::detail::coordinate_matrix. The iterator returns a transformed
 * seqan3::views::iota_simd view over the row indices. In the transformation every row index is converted to a
 * seqan3::detail::matrix_index using the seqan3::detail::convert_to_matrix_coordinate function object.
 */
template <simd_concept index_t>
class coordinate_matrix_simd_iterator
{
private:
    //!\brief The scalar type.
    using scalar_type = typename simd_traits<index_t>::scalar_type;

    //!\brief The column index as simd vector.
    index_t simd_column_index{};
    //!\brief The currently represented column index.
    scalar_type scalar_column_index{};
    //!\brief The end index of the row.
    scalar_type column_size{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(views::iota_simd<index_t>{static_cast<scalar_type>(0), column_size}
                              | std::views::transform(convert_to_matrix_coordinate<index_t>{}));
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
    coordinate_matrix_simd_iterator() noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator(coordinate_matrix_simd_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator(coordinate_matrix_simd_iterator &&) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator & operator=(coordinate_matrix_simd_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator & operator=(coordinate_matrix_simd_iterator &&) noexcept = default; //!< Defaulted.
    ~coordinate_matrix_simd_iterator() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the current column index and the row index marking the end
     *        of the rows (size of one column).
     *
     * \param[in] column_index The column index to point to.
     * \param[in] column_size The size of the column.
     */
    explicit coordinate_matrix_simd_iterator(scalar_type const column_index, scalar_type const column_size) noexcept :
        simd_column_index{simd::fill<index_t>(column_index)},
        scalar_column_index{column_index},
        column_size{column_size}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    reference operator*() const
    {
        return views::iota_simd<index_t>{static_cast<scalar_type>(0), column_size}
             | std::views::transform(convert_to_matrix_coordinate<index_t>{simd_column_index});
    }
    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    coordinate_matrix_simd_iterator & operator++()
    {
        ++scalar_column_index;
        ++simd_column_index;
        return *this;
    }

    //!\brief Increments the iterator to the next column and returns the iterator pointing to the previous one.
    coordinate_matrix_simd_iterator operator++(int)
    {
        coordinate_matrix_simd_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(coordinate_matrix_simd_iterator const & lhs, coordinate_matrix_simd_iterator const & rhs)
    {
        return lhs.scalar_column_index == rhs.scalar_column_index;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(coordinate_matrix_simd_iterator const & lhs, coordinate_matrix_simd_iterator const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}
};

} // namespace seqan3::detail
