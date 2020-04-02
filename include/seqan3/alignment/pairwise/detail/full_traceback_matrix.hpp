// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::full_traceback_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/std/span>

namespace seqan3::detail
{

template <typename trace_t>
class full_traceback_matrix
{
private:
    using vertical_column_t = decltype(views::repeat_n(trace_t{}, 1));

    class iterator;

    std::vector<trace_t> optimal_matrix{};
    std::vector<trace_t> horizontal_column{};
    vertical_column_t vertical_column{};
    size_t columns_count{};
    size_t rows_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    full_traceback_matrix() = default; //!< Defaulted.
    full_traceback_matrix(full_traceback_matrix const &) = default; //!< Defaulted.
    full_traceback_matrix(full_traceback_matrix &&) noexcept = default; //!< Defaulted.
    full_traceback_matrix & operator=(full_traceback_matrix const &) = default; //!< Defaulted.
    full_traceback_matrix & operator=(full_traceback_matrix &&) noexcept = default; //!< Defaulted.
    ~full_traceback_matrix() noexcept = default; //!< Defaulted.
    //!\}

    template <std::unsigned_integral column_count_t, std::unsigned_integral row_count_t>
    void reset_matrix(column_count_t const column_count, row_count_t const row_count)
    {
        optimal_matrix.clear();
        optimal_matrix.resize(column_count * row_count, trace_directions::none);
        horizontal_column.resize(row_count, trace_directions::none);
        vertical_column = views::repeat_n(trace_directions::none, row_count);
        this->columns_count = column_count;
        this->rows_count = row_count;
    }

    iterator begin()
    {
        return iterator{*this};
    }

    std::ranges::default_sentinel_t end()
    {
        return std::ranges::default_sentinel;
    }
};

template <typename trace_t>
class full_traceback_matrix<trace_t>::iterator
{
private:
    using combined_column_t = decltype(views::zip(std::declval<std::span<trace_t>>(),
                                                  std::declval<std::vector<trace_t> &>(),
                                                  std::declval<vertical_column_t &>()));

    static constexpr auto to_cell_proxy = std::views::transform([] (auto && tpl)
    {
        using fwd_tuple_t = decltype(tpl);
        return affine_trace_proxy<remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
    });

    full_traceback_matrix * host_ptr{nullptr};
    size_t current_column{};
    size_t end_column{};
public:

    /*!\name Associated types
        * \{
        */
    using value_type = decltype(std::declval<combined_column_t>() | to_cell_proxy);
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

    explicit iterator(full_traceback_matrix & matrix_obj) noexcept :
        host_ptr{std::addressof(matrix_obj)},
        end_column{host_ptr->columns_count}
    {}
    //!\}

    auto operator*() const
    {
        auto column_offset = std::ranges::next(host_ptr->optimal_matrix.begin(), current_column * host_ptr->rows_count);
        auto optimum_column = std::span{column_offset, std::ranges::next(column_offset, host_ptr->rows_count)};

        return views::zip(optimum_column, host_ptr->horizontal_column, host_ptr->vertical_column) | to_cell_proxy;
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

static_assert(std::ranges::sized_range<std::iter_reference_t<std::ranges::iterator_t<full_traceback_matrix<trace_directions>>>>, "Not a sized range?");
} // namespace seqan3::detail
