// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_coordinate_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// Make borrowed range
template <bool coordinate_enabled>
class alignment_coordinate_matrix
{
private:
    size_t column_count{};
    size_t row_count{};
public:

    struct alignment_coordinate : public std::pair<size_t, size_t>
    {
        using base_t = std::pair<size_t, size_t>;
        using base_t::base_t;
        using base_t::operator=;

        size_t column_index() const noexcept
        {
            return first;
        }

        size_t row_index() const noexcept
        {
            return second;
        }
    };

    class iterator;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_coordinate_matrix() = default; //!< Defaulted.
    alignment_coordinate_matrix(alignment_coordinate_matrix const &) = default; //!< Defaulted.
    alignment_coordinate_matrix(alignment_coordinate_matrix &&) = default; //!< Defaulted.
    alignment_coordinate_matrix & operator=(alignment_coordinate_matrix const &) = default; //!< Defaulted.
    alignment_coordinate_matrix & operator=(alignment_coordinate_matrix &&) = default; //!< Defaulted.
    ~alignment_coordinate_matrix() = default; //!< Defaulted.
    //!\}

    template <std::unsigned_integral column_count_t, std::unsigned_integral row_count_t>
    void reset_matrix(column_count_t column_count, row_count_t row_count)
    {
        this->column_count = column_count;
        this->row_count = row_count;
    }

    iterator begin() const noexcept
    {
        return iterator{*this};
    }

    std::ranges::default_sentinel_t end() const noexcept
    {
        return std::ranges::default_sentinel;
    }
};

template <bool coordinate_enabled>
class alignment_coordinate_matrix<coordinate_enabled>::iterator
{
private:
    struct to_alignment_coordinate
    {
        size_t m_current_column{0};

        auto operator()(size_t const row_index) noexcept
        {
            return alignment_coordinate{m_current_column, row_index};
        }
    };

    size_t current_column{0};
    size_t end_column{0};
    size_t end_row{0};

public:

    /*!\name Associated types
     * \{
     */
    using value_type = std::conditional_t<coordinate_enabled,
                                          decltype(std::views::iota(0u, end_row)
                                                 | std::views::transform(to_alignment_coordinate{})),
                                          decltype(views::repeat_n(std::ignore, 0))>;
    using reference = value_type;
    using pointer = void;
    using difference_type = std::ptrdiff_t;
    using iterator_tag = std::forward_iterator_tag;
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

    explicit iterator(alignment_coordinate_matrix const & matrix_obj) noexcept :
        end_column{matrix_obj.column_count},
        end_row{matrix_obj.row_count}
    {}
    //!\}

    reference operator*() const
    {
        if constexpr (coordinate_enabled)
            return std::views::iota(0u, end_row) | std::views::transform(to_alignment_coordinate{current_column});
        else
            return views::repeat_n(std::ignore, end_row);
    }

    iterator & operator++() noexcept
    {
        ++current_column;
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
        return lhs.current_column == rhs.current_column;
    }

    friend bool operator!=(iterator const & lhs, iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator==(iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        return lhs.current_column == lhs.end_column;
    }

    friend bool operator==(std::ranges::default_sentinel_t const & lhs, iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    friend bool operator!=(iterator const & lhs, std::ranges::default_sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    friend bool operator!=(std::ranges::default_sentinel_t const & lhs, iterator const & rhs) noexcept
    {
        return rhs != lhs;
    }
};

} // namespace seqan3::detail
