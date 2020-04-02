// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::combined_alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
template <typename ...matrices_t>
class combined_alignment_matrix : public std::ranges::view_interface<combined_alignment_matrix<matrices_t...>>
{
private:

    using combined_matrix_t = decltype(views::zip(std::declval<matrices_t>()...));

    class iterator;

    combined_matrix_t m_combined_matrix{};

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    combined_alignment_matrix() = default; //!< Defaulted.
    combined_alignment_matrix(combined_alignment_matrix const &) = default; //!< Defaulted.
    combined_alignment_matrix(combined_alignment_matrix &&) noexcept(noexcept(std::is_nothrow_move_constructible_v<combined_matrix_t>)) = default; //!< Defaulted.
    combined_alignment_matrix & operator=(combined_alignment_matrix const &) = default; //!< Defaulted.
    combined_alignment_matrix & operator=(combined_alignment_matrix &&) noexcept(noexcept(std::is_nothrow_move_assignable_v<combined_matrix_t>)) = default; //!< Defaulted.
    ~combined_alignment_matrix() noexcept = default; //!< Defaulted.

    combined_alignment_matrix(matrices_t ...matrices) :
        m_combined_matrix{views::zip(std::forward<matrices_t>(matrices)...)}
    {}
    //!\}

    iterator begin()
    {
        return iterator{m_combined_matrix};
    }

    std::ranges::sentinel_t<combined_matrix_t> end()
    {
        return std::ranges::end(m_combined_matrix);
    }
};

template <std::ranges::viewable_range ... matrices_t>
class combined_alignment_matrix<matrices_t...>::iterator
{
private:
    using combined_matrix_iterator_type = std::ranges::iterator_t<combined_matrix_t>;
    using combined_matrix_sentinel_type = std::ranges::sentinel_t<combined_matrix_t>;
    using combined_column_type = decltype(views::zip(std::declval<std::ranges::range_reference_t<matrices_t>>()...));

    combined_matrix_iterator_type m_combined_matrix_iterator{};
public:
    /*!\name Associated types
     * \{
     */
    using value_type = std::conditional_t<sizeof...(matrices_t) == 1,
                                         std::tuple_element_t<0, std::ranges::range_value_t<combined_matrix_t>>,
                                         combined_column_type>;
    using reference = std::conditional_t<sizeof...(matrices_t) == 1,
                                         std::tuple_element_t<0, std::ranges::range_reference_t<combined_matrix_t>>,
                                         combined_column_type>;
    using pointer = void;
    using difference_type = std::ptrdiff_t;
    using iterator_tag = std::forward_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    iterator() = default; //!< Defaulted.
    iterator(iterator const &) = default; //!< Defaulted.
    iterator(iterator &&) = default; //!< Defaulted.
    iterator & operator=(iterator const &) = default; //!< Defaulted.
    iterator & operator=(iterator &&) = default; //!< Defaulted.
    ~iterator() noexcept = default; //!< Defaulted.

    iterator(combined_matrix_t & combined_matrix) noexcept :
        m_combined_matrix_iterator{std::ranges::begin(combined_matrix)}
    {}
    //!\}

    reference operator*() const
    {
        return std::apply([] (auto && ...columns)
        {
            if constexpr (sizeof...(columns) == 1)
                return (std::forward<decltype(columns)>(columns), ...);
            else
                return views::zip(std::forward<decltype(columns)>(columns)...);
        }, *m_combined_matrix_iterator);
    }

    iterator & operator++() noexcept
    {
        ++m_combined_matrix_iterator;
        return *this;
    }

    iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    friend bool operator==(iterator const & lhs, iterator const & rhs) noexcept
    {
        return lhs.m_combined_matrix_iterator == rhs.m_combined_matrix_iterator;
    }

    friend bool operator!=(iterator const & lhs, iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator==(iterator const & lhs, combined_matrix_sentinel_type const & rhs) noexcept
    {
        return lhs.m_combined_matrix_iterator == rhs;
    }

    friend bool operator==(combined_matrix_sentinel_type const & lhs, iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    friend bool operator!=(iterator const & lhs, combined_matrix_sentinel_type const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(combined_matrix_sentinel_type const & lhs, iterator const & rhs) noexcept
    {
        return rhs != lhs;
    }
};
} // namespace seqan3::detail
