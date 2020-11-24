// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_storage_two_dimensional.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_concept.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>

namespace seqan3::detail
{

// Concept:


template <alignment_matrix_element matrix_element_t>
class alignment_matrix_storage_two_dimensional : protected matrix_element_t
{
private:
    //!\brief The combined column score.
    using storage_value_type = typename matrix_element_t::storage_value_type;
    //!\brief The type of the score column which allocates memory for the entire column.
    using storage_type = two_dimensional_matrix<storage_value_type,
                                                aligned_allocator<storage_value_type, alignof(storage_value_type)>,
                                                matrix_major_order::column>;

    //!\brief The column over the optimal scores.
    storage_type matrix{};

    class iterator;

public:

    using matrix_element_type = matrix_element_t;

    // We might be able to optimise this to a span!
    // It depends on the element_type.
    using column_type = std::conditional_t<matrix_element_t::returns_element_proxy,
                                           std::ranges::subrange<iterator, iterator>,
                                           std::span<storage_value_type>>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_matrix_storage_two_dimensional() = default; //!< Defaulted.
    alignment_matrix_storage_two_dimensional(alignment_matrix_storage_two_dimensional const &)
        = default; //!< Defaulted.
    alignment_matrix_storage_two_dimensional(alignment_matrix_storage_two_dimensional &&) = default; //!< Defaulted.
    alignment_matrix_storage_two_dimensional & operator=(alignment_matrix_storage_two_dimensional const &)
        = default; //!< Defaulted.
    alignment_matrix_storage_two_dimensional & operator=(alignment_matrix_storage_two_dimensional &&)
        = default; //!< Defaulted.
    ~alignment_matrix_storage_two_dimensional() = default; //!< Defaulted.
    //!\}

    template <std::integral column_index_t, std::integral row_index_t>
    void resize(row_index_type<row_index_t> const row_count,
                column_index_type<column_index_t> const column_count,
                typename matrix_element_t::value_type const initial_value = typename matrix_element_t::value_type{})
    {
        // strong exception guarantee!
        matrix.resize(number_rows{row_count.get()}, number_cols{column_count.get()}, this->initialise(initial_value));
    }

    void clear()
    {
        matrix.clear();
    }

    size_t cols() const noexcept
    {
        return matrix.cols();
    }

    size_t rows() const noexcept
    {
        return matrix.rows();
    }

    column_type column_at(size_t const column_id)
    {
        auto it = matrix.begin() + matrix_offset{row_index_type{0}, column_index_type{static_cast<int64_t>(column_id)}};
        if constexpr (matrix_element_t::returns_element_proxy)
            return column_type{iterator{this, it}, iterator{this, std::ranges::next(it, rows())}};
        else
            return column_type{std::addressof(*it), rows()};
    }
};

// We only need the inherited iterator and overload the dereference operator.
template <typename matrix_element_t>
class alignment_matrix_storage_two_dimensional<matrix_element_t>::iterator :
    public inherited_iterator_base<iterator, std::ranges::iterator_t<storage_type>>
{
private:
    using underlying_iterator_type = std::ranges::iterator_t<storage_type>;
    using base_type = inherited_iterator_base<iterator, underlying_iterator_type>;

    friend base_type;

    alignment_matrix_storage_two_dimensional * host_ptr{nullptr};
public:
    /*!\name Associated types
     * \{
     */
    // typename concept::element_type
    //!\brief The value type.
    using value_type = typename matrix_element_t::element_type; //element<std::tuple<score_t, score_t, score_t>>;
    // typename concept::element_reference
    //!\brief The reference type.
    using reference = typename matrix_element_t::element_reference; // element<std::tuple<score_t &, score_t &, score_t &>>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using typename base_type::difference_type;
    //!\brief The iterator category.
    using typename base_type::iterator_category;
    //!\brief The iterator concept.
    using typename base_type::iterator_concept;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    iterator() = default; //!< Defaulted.
    iterator(iterator const &) = default; //!< Defaulted.
    iterator(iterator &&) = default; //!< Defaulted.
    iterator & operator=(iterator const &) = default; //!< Defaulted.
    iterator & operator=(iterator &&) = default; //!< Defaulted.
    ~iterator() = default; //!< Defaulted.

    iterator(alignment_matrix_storage_two_dimensional * host_ptr, underlying_iterator_type underlying_it) noexcept :
        base_type{std::move(underlying_it)},
        host_ptr{host_ptr}
    {}
    //!\}

    /*!\brief Access operators
     * \{
     */
    reference operator*() const noexcept
    {
        return host_ptr->make_element(*static_cast<base_type const &>(*this));
    }

    reference operator->() const noexcept = delete;

    reference operator[](difference_type const skip) const noexcept
    {
        return *(*this + skip);
    }
    //!\}

private:

    iterator from_base_impl(underlying_iterator_type base_it) const
    {
        return iterator{host_ptr, base_it};
    }
};

}  // namespace seqan3::detail

// So if we only consider the pure score matrix we do not need to do anything.
// We do not want to rewrite the entire matrix just because we need an additional value
// We do not want to pay for another value.
// So we add a level of indirection
// This way we can control the vertical value.
// We need to manage the element type though.
// The element type is than presented by the concept.

// So the alignment matrix defines a column based iteration over a matrix
// It is created with a concept for the alignment_matrix_storage
// alignment_matrix<alignment_matrix_storage[, alignment_matrix_storage]> -> implements alignment matrix concept
    // column_at(size_t) -> returns column
    // constructed from alignment_matrix_storage

    // alignment_matrix_storage_two_dimensional<matrix_element> -> implements alignment_matrix_storage concept
        // this tells me that the storage is based on a single column
            // we do optimise it by also allowing to store only one vertical cell for the score
            // in another case we need to make the vertical cell hold two values.
        // score_affine_light -> implements matrix_element concept
        // score_affine_trace -> implements matrix_element concept
        // score_linear -> implements matrix_element concept
        //
    // alignment_matrix_storage_two_dimensional<matrix_element>
        // matrix_element_single -> single value can be a trace.
    // alignment_matrix_storage_sampling<matrix_element>


// alignment matrix only implements the column idea:
// must be able to ask the number of columns of the underlying storage
