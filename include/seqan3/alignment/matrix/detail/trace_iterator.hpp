// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_base.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A CRTP-base class for trace iterator implementations for the alignment algorithms.
 * \ingroup alignment_matrix
 * \implements std::forward_iterator
 *
 * \tparam matrix_iter_t The wrapped matrix iterator; must model seqan3::detail::two_dimensional_matrix_iterator and
 *                       the iterator's value type must be the same as seqan3::detail::trace_directions, i.e.
 *                       `std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>` must evaluate to `true`.
 *
 * \details
 *
 * This iterator follows the trace path generated by the alignment algorithm. It wraps an underlying
 * seqan3::detail::two_dimensional_matrix_iterator over a trace matrix, whose value type is
 * seqan3::detail::trace_directions. The iterator moves along the trace path until it finds a cell with
 * seqan3::detail::trace_directions::none.
 * Accordingly, when advancing this iterator, it actually moves from right to left and from bottom to top in the
 * underlying matrix. When the iterator is dereferenced, it outputs any of the following direction:
 * seqan3::detail::trace_directions::diagonal, seqan3::detail::trace_directions::up, or
 * seqan3::detail::trace_directions::left.
 *
 * In addition, the iterator provides an additional member to access the current position as a
 * seqan3::detail::matrix_coordinate.
 *
 * This iterator also models the [Cpp17ForwardIterator](https://en.cppreference.com/w/cpp/named_req/ForwardIterator).
 * Note, it does not directly dereference the actual trace direction stored in the underlying matrix.
 * Thus, it cannot be used as an output iterator.
 *
 * ### Overloading the behaviour
 *
 * The behaviour of following a trace direction can be customised through the derived type by overloading the functions
 * * seqan3::detail::trace_iterator::go_diagonal,
 * * seqan3::detail::trace_iterator::go_left, and
 * * seqan3::detail::trace_iterator::go_up.
 *
 * In the default implementation they move along an unbanded matrix. This means, they go to the previous cell in the
 * respective direction.
 */
template <two_dimensional_matrix_iterator matrix_iter_t>
class trace_iterator
{
private:
    static_assert(std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>,
                  "Value type of the underlying iterator must be seqan3::detail::trace_directions.");

    template <two_dimensional_matrix_iterator>
    friend class trace_iterator;

    matrix_iter_t matrix_iter{}; //!< The underlying matrix iterator.
    size_t pivot_column{}; //!< The largest column index which is inside of the band in the first row of the matrix.
    trace_directions current_direction{}; //!< The current trace direction.
    bool legacy_iterator{false}; //!< If this iterator is used by the old alignment algorithm.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr trace_iterator() = default; //!< Defaulted.
    constexpr trace_iterator(trace_iterator const &) = default; //!< Defaulted.
    constexpr trace_iterator(trace_iterator &&) = default; //!< Defaulted.
    constexpr trace_iterator & operator=(trace_iterator const &) = default; //!< Defaulted.
    constexpr trace_iterator & operator=(trace_iterator &&) = default; //!< Defaulted.
    ~trace_iterator() = default; //!< Defaulted.

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \param[in] matrix_iter The underlying matrix iterator.
     */
    explicit constexpr trace_iterator(matrix_iter_t matrix_iter) noexcept :
        trace_iterator{std::move(matrix_iter), column_index_type{std::numeric_limits<size_t>::max()}, false}
    {}

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \param[in] matrix_iter The underlying matrix iterator.
     * \param[in] pivot_column The last column index which is still inside of the band in the first row of the
     *                         banded matrix.
     * \param[in] legacy_iterator A boolean flag indicating whether the old alignment is using this; defaults to `true`.
     *
     * \details
     *
     * The legacy_iterator flag changes the behaviour of the banded iterator when moving to the previous cell in
     * horizontal or diagonal direction. This depends on the alignment implementation and is handled differently in
     * the new implementation.
     */
    template <typename index_t>
    constexpr trace_iterator(matrix_iter_t matrix_iter,
                             column_index_type<index_t> const & pivot_column,
                             bool legacy_iterator = true)
        noexcept :
        matrix_iter{std::move(matrix_iter)},
        pivot_column{static_cast<size_t>(pivot_column.get())},
        legacy_iterator{legacy_iterator}
    {
        set_trace_direction(*this->matrix_iter);
    }

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \tparam other_matrix_iter_t The underlying matrix iterator type of `other`; the condition
     *                             `std::constructible_from<matrix_iter_t, other_matrix_iter_t>` must evaluate to `true`.
     * \param[in] other The underlying matrix iterator.
     *
     * \details
     *
     * Allows the conversion of non-const to const iterator.
     */
    template <two_dimensional_matrix_iterator other_matrix_iter_t>
    //!\cond
        requires (!std::same_as<matrix_iter_t, other_matrix_iter_t>) &&
                 std::constructible_from<matrix_iter_t, other_matrix_iter_t>
    //!\endcond
    constexpr trace_iterator(trace_iterator<other_matrix_iter_t> other) noexcept :
        trace_iterator{matrix_iter_t{std::move(other.matrix_iter)},
                       column_index_type{other.pivot_column},
                       other.legacy_iterator}
    {}

    //!\}

    /*!\name Associated types
     * \{
     */
    using value_type = trace_directions; //!< The value type.
    using reference = trace_directions const &; //!< The reference type.
    using pointer = value_type *; //!< The pointer type.
    using difference_type = std::ptrdiff_t; //!< The difference type.
    using iterator_category = std::forward_iterator_tag; //!< Forward iterator tag.
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the current trace direction.
    reference operator*() const noexcept
    {
        return current_direction;
    }

    //!\brief Returns a pointer to the current trace direction.
    pointer operator->() const noexcept
    {
        return &current_direction;
    }

    //!\brief Returns the current coordinate in two-dimensional space.
    [[nodiscard]] constexpr matrix_coordinate coordinate() const noexcept
    {
        auto coord = matrix_iter.coordinate();
        // only correct the row coordinate if the current column is greater than the set upper_diagonal.
        coord.row += (coord.col > pivot_column || legacy_iterator) * static_cast<int32_t>(coord.col - pivot_column);
        return coord;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advances the iterator by one.
    constexpr trace_iterator & operator++() noexcept
    {
        trace_directions old_dir = *matrix_iter;

        assert(old_dir != trace_directions::none);

        if (current_direction == trace_directions::up)
        {
            go_up(matrix_iter);
            // Set new trace direction if last position was up_open.
            if (static_cast<bool>(old_dir & trace_directions::up_open))
                set_trace_direction(*matrix_iter);
        }
        else if (current_direction == trace_directions::left)
        {
            go_left(matrix_iter);
            // Set new trace direction if last position was left_open.
            if (static_cast<bool>(old_dir & trace_directions::left_open))
                set_trace_direction(*matrix_iter);
        }
        else
        {
            assert(current_direction == trace_directions::diagonal);

            go_diagonal(matrix_iter);
            set_trace_direction(*matrix_iter);
        }
        return *this;
    }

    //!\brief Returns an iterator advanced by one.
    constexpr trace_iterator operator++(int) noexcept
    {
        trace_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Returns `true` if both iterators are equal, `false` otherwise.
    constexpr friend bool operator==(trace_iterator const & lhs, trace_iterator const & rhs) noexcept
    {
        return lhs.matrix_iter == rhs.matrix_iter;
    }

    //!\brief Returns `true` if the pointed-to-element is seqan3::detail::trace_directions::none.
    constexpr friend bool operator==(trace_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return *lhs.matrix_iter == trace_directions::none;
    }

    //!\brief copydoc operator==()
    constexpr friend bool operator==(std::default_sentinel_t const &, trace_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief Returns `true` if both iterators are not equal, `false` otherwise.
    constexpr friend bool operator!=(trace_iterator const & lhs, trace_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Returns `true` if the pointed-to-element is not seqan3::detail::trace_directions::none.
    constexpr friend bool operator!=(trace_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief copydoc operator!=()
    constexpr friend bool operator!=(std::default_sentinel_t const &, trace_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}

private:
    /*!\name Overload functions
     * \brief These functions can be overloaded by the derived class to customise the iterator.
     * \{
     */
    //!\brief Moves iterator to previous left cell.
    constexpr void go_left(matrix_iter_t & iter) const noexcept
    {
        // Note, in the banded matrix, the columns are virtually shifted by one cell.
        // So going left means go to the previous column and then one row down.
        int32_t row = coordinate().col > pivot_column || legacy_iterator;
        iter -= matrix_offset{row_index_type{-row}, column_index_type{1}};
    }

    //!\brief Moves iterator to previous up cell.
    constexpr void go_up(matrix_iter_t & iter) const noexcept
    {
        iter -= matrix_offset{row_index_type{1}, column_index_type{0}};
    }

    //!\brief Moves iterator to previous diagonal cell.
    constexpr void go_diagonal(matrix_iter_t & iter) const noexcept
    {
        // Note, in the banded matrix, the columns are virtually shifted by one cell.
        // So going diagonal means go to the previous column and stay in the same row.
        int32_t row = coordinate().col <= pivot_column && !legacy_iterator;
        iter -= matrix_offset{row_index_type{row}, column_index_type{1}};
    }
    //!\}

    //!\brief Updates the current trace direction.
    void set_trace_direction(trace_directions const dir) noexcept
    {
        if (static_cast<bool>(dir & trace_directions::diagonal))
        {
            current_direction = trace_directions::diagonal;
        }
        else if (static_cast<bool>(dir & trace_directions::up) ||
                 static_cast<bool>(dir & trace_directions::up_open))
        {
            current_direction = trace_directions::up;
        }
        else if (static_cast<bool>(dir & trace_directions::left) ||
                 static_cast<bool>(dir & trace_directions::left_open))
        {
            current_direction = trace_directions::left;
        }
        else
        {
            current_direction = trace_directions::none;
        }
    }
};

} // namespace seqan3::detail
