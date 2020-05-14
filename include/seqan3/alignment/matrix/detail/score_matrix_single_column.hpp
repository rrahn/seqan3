// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::score_matrix_single_column.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Score matrix for the pairwise alignment using only a single column.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_t The type of the score; must model seqan3::arithmetic.
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
template <typename score_t>
//!\cond
    requires arithmetic<score_t> || simd_concept<score_t>
//!\endcond
class score_matrix_single_column
{
private:
    //!\brief The type of the score column which allocates memory for the entire column.
    using physical_column_t = std::vector<score_t, aligned_allocator<score_t, alignof(score_t)>>;
    //!\brief The type of the virtual score column which only stores one value.
    using virtual_column_t = decltype(views::repeat_n(score_t{}, 1));

    class matrix_iterator;

    //!\brief The column over the optimal scores.
    physical_column_t m_optimal_column{};
    //!\brief The column over the horizontal gap scores.
    physical_column_t m_horizontal_column{};
    //!\brief The virtual column over the vertical gap scores.
    virtual_column_t m_vertical_column{};
    //!\brief The number of columns for this matrix.
    size_t m_number_of_columns{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    score_matrix_single_column() = default; //!< Defaulted.
    score_matrix_single_column(score_matrix_single_column const &) = default; //!< Defaulted.
    score_matrix_single_column(score_matrix_single_column &&) = default; //!< Defaulted.
    score_matrix_single_column & operator=(score_matrix_single_column const &) = default; //!< Defaulted.
    score_matrix_single_column & operator=(score_matrix_single_column &&) = default; //!< Defaulted.
    ~score_matrix_single_column() = default; //!< Defaulted.
    //!\}

    /*!\brief Resets the matrix dimensions given by two sequences.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     *
     * \details
     *
     * Resizes the optimal and the horizontal score column to the new dimensions given by the sizes of the
     * two sequences. The number of columns is set to the size of sequence1 plus 1 for the initialisation and
     * the number of rows is set to the size of sequence2 plus 1 for the initialisation as well.
     * Reallocation happens only if the new column size exceeds the current capacity of the optimal and horizontal
     * score column. The underlying vectors will not be cleared during the reset.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void reset_matrix(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        m_number_of_columns = std::ranges::distance(sequence1) + 1;
        size_t number_of_rows = std::ranges::distance(sequence2) + 1;
        m_optimal_column.resize(number_of_rows);
        m_horizontal_column.resize(number_of_rows);
        m_vertical_column = views::repeat_n(score_t{}, number_of_rows);
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column.
    matrix_iterator begin()
    {
        return matrix_iterator{*this};
    }
    //!\brief This score matrix is not const-iterable.
    matrix_iterator begin() const = delete;

    //!\brief Returns a default sentinel indicating the end of the matrix.
    std::ranges::default_sentinel_t end()
    {
        return std::ranges::default_sentinel;
    }

    //!\brief This score matrix is not const-iterable.
    std::ranges::default_sentinel_t end() const = delete;
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
template <typename score_t>
class score_matrix_single_column<score_t>::matrix_iterator
{
private:

    //!\brief The type of the zipped score column.
    using matrix_column_t = decltype(views::zip(std::declval<physical_column_t &>(),
                                                std::declval<physical_column_t &>(),
                                                std::declval<virtual_column_t &>()));

    //!\brief The transform adaptor to convert the tuple from the zip view into a seqan3::detail::affine_cell_type.
    static constexpr auto transform_to_affine_cell = std::views::transform([] (auto && tpl)
    {
        using fwd_tuple_t = decltype(tpl);
        return affine_cell_proxy<remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
    });

    //!\brief The pointer to the underlying matrix.
    score_matrix_single_column * m_host_ptr{nullptr};
    //!\brief The current column index.
    size_t m_current_column{};
    //!\brief The column index of the column behind the last one.
    size_t m_end_column{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::declval<matrix_column_t>() | transform_to_affine_cell);
    //!\brief The reference type.
    using reference = value_type;
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
    matrix_iterator() noexcept = default; //!< Defaulted.
    matrix_iterator(matrix_iterator const &) noexcept = default; //!< Defaulted.
    matrix_iterator(matrix_iterator &&) noexcept = default; //!< Defaulted.
    matrix_iterator & operator=(matrix_iterator const &) noexcept = default; //!< Defaulted.
    matrix_iterator & operator=(matrix_iterator &&) noexcept = default; //!< Defaulted.
    ~matrix_iterator() = default; //!< Defaulted.

    /*!\brief Initialises the iterator from the underlying matrix.
     *
     * \param[in] host_matrix The underlying matrix.
     */
    explicit matrix_iterator(score_matrix_single_column & host_matrix) noexcept :
        m_host_ptr{std::addressof(host_matrix)},
        m_end_column{m_host_ptr->m_number_of_columns}
    {}
    //!\}

    /*!\name Access operators
     * \{
     */
    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        return views::zip(m_host_ptr->m_optimal_column, m_host_ptr->m_horizontal_column, m_host_ptr->m_vertical_column)
             | transform_to_affine_cell;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    matrix_iterator & operator++()
    {
        ++m_current_column;
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
    //!\brief Compares the matrix iterator with the sentinel.
    friend bool operator==(matrix_iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        return lhs.m_current_column == lhs.m_end_column;
    }

    //!\brief Compares the matrix iterator with the sentinel.
    friend bool operator==(std::ranges::default_sentinel_t const & lhs, matrix_iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\brief Compares the matrix iterator with the sentinel.
    friend bool operator!=(matrix_iterator const & lhs, std::ranges::default_sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compares the matrix iterator with the sentinel.
    friend bool operator!=(std::ranges::default_sentinel_t const & lhs, matrix_iterator const & rhs) noexcept
    {
        return rhs != lhs;
    }
    //!\}
};

} // namespace seqan3::detail
