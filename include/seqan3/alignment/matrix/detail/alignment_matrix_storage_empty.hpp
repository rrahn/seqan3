// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_storage_empty.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>

namespace seqan3::detail
{

template <typename matrix_element_t>
class alignment_matrix_storage_empty : protected matrix_element_t
{
private:
    //!\brief The combined column score.
    using storage_value_type = typename matrix_element_t::storage_value_type;

protected:

    using column_type = std::views::empty<storage_value_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_matrix_storage_empty() = default; //!< Defaulted.
    alignment_matrix_storage_empty(alignment_matrix_storage_empty const &) = default; //!< Defaulted.
    alignment_matrix_storage_empty(alignment_matrix_storage_empty &&) = default; //!< Defaulted.
    alignment_matrix_storage_empty & operator=(alignment_matrix_storage_empty const &) = default; //!< Defaulted.
    alignment_matrix_storage_empty & operator=(alignment_matrix_storage_empty &&) = default; //!< Defaulted.
    ~alignment_matrix_storage_empty() = default; //!< Defaulted.
    //!\}

    template <std::integral column_index_t, std::integral row_index_t>
    void resize(row_index_type<row_index_t> const /*row_count*/,
                column_index_type<column_index_t> const /*column_count*/,
                typename matrix_element_t::value_type const /*initial_value*/ = typename matrix_element_t::value_type{})
    {}

    void clear()
    {}

    size_t cols() const noexcept
    {
        return 0;
    }

    size_t rows() const noexcept
    {
        return 0;
    }

    column_type column_at(size_t const /*column_id*/)
    {
        return column_type{};
    }
};
