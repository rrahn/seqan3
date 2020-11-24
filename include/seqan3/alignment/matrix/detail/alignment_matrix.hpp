// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

/*!\brief Score matrix for the pairwise alignment using only a single column.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_t The type of the score; must model seqan3::arithmetic or seqan3::simd::simd_concept.
 *
 * \details
 *
 * In many cases it is sufficient to store only a single score column to compute the alignment between two sequences.
 * Since the alignment is computed iteratively column by column, the same memory can be reused for the next score.
 * This score matrix stores the complete column for both the optimal and horizontal score, but only stores a single
 * value for the vertical column. Hence, this matrix can only be used for a column
 * based computation layout.
 *
 * ### Range interface
 *
 * The matrix offers a input range interface over the columns of the matrix. Dereferencing the iterator will return
 * another range which represents the actual score column in memory. The returned range is a
 * transformed seqan3::views::zip view over the optimal, horizontal and vertical column. The reference type of this
 * view is the seqan3::detail::affine_cell_proxy, which offers a practical interface to access the value of the
 * optimal, horizontal and vertical value of the underlying matrices.
 */
template <typename matrix_storage_t>
class alignment_matrix
{
private:

    template <typename matrix_element_t>
    class column_value;

    class iterator;

    using column_type = decltype(std::declval<matrix_storage_t &>().column_at(0u));

    matrix_storage_t matrix_storage{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_matrix() = default; //!< Defaulted.
    alignment_matrix(alignment_matrix const &) = default; //!< Defaulted.
    alignment_matrix(alignment_matrix &&) = default; //!< Defaulted.
    alignment_matrix & operator=(alignment_matrix const &) = default; //!< Defaulted.
    alignment_matrix & operator=(alignment_matrix &&) = default; //!< Defaulted.
    ~alignment_matrix() = default; //!< Defaulted.

    /*!\brief Constructs a new alignment matrix with the given storage.
     *
     * \param[in] storage The storage to use for this matrix.
     *
     * \details
     *
     * ### Exception
     *
     * May throw implementation-defined exceptions.
     * If an exception is thrown, this function has no effect (strong exception guarantee).
     */
    explicit alignment_matrix(matrix_storage_t storage) : matrix_storage{std::move(storage)}
    {}
    //!\}

    /*!\brief Sets the associated matrix storage to the given storage.
     *
     * \param[in] storage The storage to set.
     *
     * \details
     *
     * ### Exception
     *
     * May throw implementation-defined exceptions.
     * If an exception is thrown, this function has no effect (strong exception guarantee).
     */
    void storage(matrix_storage_t storage)
    {
        matrix_storage = std::move(storage);
    }

    /*!\brief Returns the associated matrix storage.
     *
     * \details
     *
     * Returns a reference to the associated matrix storage.
     *
     * ### Exception
     *
     * Nothrow exception guarantee.
     */
    matrix_storage_t & storage() noexcept
    {
        return matrix_storage;
    }

    //!\overload
    matrix_storage_t const & storage() const noexcept
    {
        return matrix_storage;
    }

    size_t cols() const noexcept
    {
        return matrix_storage.cols();
    }

    size_t rows() const noexcept
    {
        return matrix_storage.rows();
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column.
    iterator begin()
    {
        return iterator{this, 0u};
    }

    //!\brief This score matrix is not const-iterable.
    iterator begin() const = delete;

    //!\brief Returns the iterator pointing behind the last column.
    iterator end()
    {
        return iterator{this, cols()};
    }

    //!\brief This score matrix is not const-iterable.
    iterator end() const = delete;
    //!\}
};

/*!\brief Score matrix iterator for the pairwise alignment using only a single column.
 * \implements std::input_iterator
 *
 * \details
 *
 * Implements a counted iterator to simulate the iteration over the actual matrix. When dereferenced, the
 * iterator returns a view over the allocated memory of the respective columns. The returned view zips
 * the three columns into a single range and transforms the returned tuple to a
 * seqan3::detail::affine_cell_proxy to simplify the access to the correct values without knowing the internal
 * tuple layout returned by the seqan3::views::zip view.
 */
template <typename matrix_storage_t>
class alignment_matrix<matrix_storage_t>::iterator
{
private:

    //!\brief The pointer to the underlying matrix.
    alignment_matrix * host_ptr{nullptr};
    //!\brief The current column index.
    size_t column_id{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The reference type.
    using reference = column_type;
    //!\brief The value type.
    using value_type = column_value<std::ranges::range_value_t<reference>>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    iterator() noexcept = default; //!< Defaulted.
    iterator(iterator const &) noexcept = default; //!< Defaulted.
    iterator(iterator &&) noexcept = default; //!< Defaulted.
    iterator & operator=(iterator const &) noexcept = default; //!< Defaulted.
    iterator & operator=(iterator &&) noexcept = default; //!< Defaulted.
    ~iterator() = default; //!< Defaulted.

    /*!\brief Initialises the iterator from the underlying matrix.
     *
     * \param[in] matrix A pointer to the underlying matrix.
     * \param[in] column_id The initial column index to set the iterator to.
     */
    explicit iterator(alignment_matrix * matrix, size_t const column_id) noexcept :
        host_ptr{matrix},
        column_id{column_id}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        return host_ptr->matrix_storage.column_at(column_id);
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    iterator & operator++()
    {
        ++column_id;
        return *this;
    }

    //!\brief Move `this` to the next column.
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(iterator const & lhs, iterator const & rhs) noexcept
    {
        return lhs.column_id == rhs.column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(iterator const & lhs, iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

template <typename matrix_storage_t>
template <typename matrix_element_t>
class alignment_matrix<matrix_storage_t>::column_value :
    public std::vector<matrix_element_t, aligned_allocator<matrix_element_t, sizeof(matrix_element_t)>>
{
private:
    using base_t = std::vector<matrix_element_t, aligned_allocator<matrix_element_t, sizeof(matrix_element_t)>>;
public:

    /*!\name Constructor, assignment and destructor
     * \{
     */
    column_value() noexcept = default; //!< Defaulted.
    column_value(column_value const &) noexcept = default; //!< Defaulted.
    column_value(column_value &&) noexcept = default; //!< Defaulted.
    column_value & operator=(column_value const &) noexcept = default; //!< Defaulted.
    column_value & operator=(column_value &&) noexcept = default; //!< Defaulted.
    ~column_value() = default; //!< Defaulted.

    column_value(column_type column) :
        base_t{std::ranges::begin(column), std::ranges::end(column)}
    {}

    column_value & operator=(column_type column)
    {
        static_cast<base_t &>(*this) = base_t{std::ranges::begin(column), std::ranges::end(column)};
        return *this;
    }
    //!\}
};

} // namespace seqan3::detail
