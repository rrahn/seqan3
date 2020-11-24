// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_storage_sparse.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_concept.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

template <alignment_matrix_element matrix_element_t>
class alignment_matrix_storage_sparse : public matrix_element_t
{
private:
    using value_type = typename matrix_element_t::value_type;
    using storage_value_type = typename matrix_element_t::storage_value_type;
    using storage_type = std::vector<storage_value_type, seqan3::aligned_allocator<storage_value_type, sizeof(storage_value_type)>>;
    using vertical_column_t = decltype(views::repeat_n(value_type{}, 0));

    class iterator;

    //!\brief The column over the optimal scores.
    storage_type sparse_matrix{};
    size_t column_count{};
public:

    using matrix_element_type = matrix_element_t;

    using column_type = std::conditional_t<matrix_element_t::returns_element_proxy,
                                           std::ranges::subrange<iterator, iterator>,
                                           std::span<storage_value_type>>;
    // using column_type = std::span<storage_value_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_matrix_storage_sparse() = default; //!< Defaulted.
    alignment_matrix_storage_sparse(alignment_matrix_storage_sparse const &) = default; //!< Defaulted.
    alignment_matrix_storage_sparse(alignment_matrix_storage_sparse &&) = default; //!< Defaulted.
    alignment_matrix_storage_sparse & operator=(alignment_matrix_storage_sparse const &) = default; //!< Defaulted.
    alignment_matrix_storage_sparse & operator=(alignment_matrix_storage_sparse &&) = default; //!< Defaulted.
    ~alignment_matrix_storage_sparse() = default; //!< Defaulted.
    //!\}

    // Here we actually require that the initial value is semiregular.
    // This must be handled in the concept definition of alignment_matrix_element
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(row_index_type<row_index_t> const row_count,
                column_index_type<column_index_t> const column_count,
                typename matrix_element_t::value_type const initial_value = typename matrix_element_t::value_type{})
    {
        // strong exception guarantee!
        sparse_matrix.resize(row_count.get(), this->initialise(initial_value));
        this->column_count = column_count.get();
    }

    void clear()
    {
        sparse_matrix.clear();
        column_count = 0;
    }

    size_t cols() const noexcept
    {
        return column_count;
    }

    size_t rows() const noexcept
    {
        return sparse_matrix.size();
    }

    auto column_at([[maybe_unused]] size_t const column_id)
    {
        assert(column_id < cols());
        if constexpr (matrix_element_t::returns_element_proxy)
            return column_type{iterator{this, sparse_matrix.begin()}, iterator{this, sparse_matrix.end()}};
        else
            return column_type{sparse_matrix};
    }
};

template <alignment_matrix_element matrix_element_t>
class alignment_matrix_storage_sparse<matrix_element_t>::iterator //:
    // public inherited_iterator_base<iterator, typename storage_type::iterator>
{
private:
    using underlying_iterator_type = std::ranges::iterator_t<storage_type>;
    // using base_type = inherited_iterator_base<iterator, underlying_iterator_type>;

    // friend base_type;

    alignment_matrix_storage_sparse * host_ptr{nullptr};
    underlying_iterator_type base_it{};
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
    using difference_type = typename underlying_iterator_type::difference_type;
    // using typename base_type::difference_type;
    //!\brief The iterator category.
    using iterator_category = std::forward_iterator_tag; //typename underlying_iterator_type::iterator_category;
    //!\brief The iterator concept.
    // using typename underlying_iterator_type::iterator_concept;
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

    iterator(alignment_matrix_storage_sparse * host_ptr, underlying_iterator_type underlying_it) noexcept :
        host_ptr{host_ptr},
        base_it{std::move(underlying_it)}
    {}
    //!\}

    /*!\brief Access operators
     * \{
     */
    reference operator*() const noexcept
    {
        assert(host_ptr != nullptr);
        return host_ptr->make_element(*base_it);
    }

    pointer operator->() const noexcept = delete;

    iterator & operator++() noexcept
    {
        ++base_it;
        return *this;
    }

    iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr bool operator==(iterator const & rhs) const noexcept
    {
        return base_it == rhs.base_it;
    }

    constexpr bool operator!=(iterator const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    // reference operator[](difference_type const skip) const noexcept
    // {
    //     return *(*this + skip);
    // }
    //!\}
private:

    /*!\brief Generates a new iterator from the base iterator.
     * \param[in] base_it The base iterator to construct from.
     */
    // iterator from_base_impl(underlying_iterator_type base_it) const noexcept
    // {
    //     return iterator{host_ptr, std::move(base_it)};
    // }
};
} // namespace seqan3::detail

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

    // alignment_matrix_storage_sparse<matrix_element> -> implements alignment_matrix_storage concept
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
