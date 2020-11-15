// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_matrix_simd_adapter_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/iterator>

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_base.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/simd/concept.hpp>

namespace seqan3::detail
{
/*!\brief A two-dimensional matrix iterator adapter that converts a simd trace value into
 *        seqan3::detail::trace_directions object for a given simd lane.
 * \ingroup alignment_matrix
 *
 * \tparam matrix_iter_t The underlying wrapped seqan3::detail::two_dimensional_matrix iterator; the iterator value type
 *                       must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Wraps a two-dimensional matrix iterator over the simd trace matrix. On dereferencing, the value at the given
 * simd lane of the underlying simd vector is convert to a seqan3::detail::trace_directions object.
 * This wrapper allows to build a trace path around a simd trace matrix and use the
 * seqan3::detail::aligned_sequence_builder without transforming the memory of the underlying trace matrix.
 *
 * \remark This iterator does not model std::output_iterator.
 */
template <two_dimensional_matrix_iterator matrix_iter_t>
//!\cond
    requires simd::simd_concept<std::iter_value_t<matrix_iter_t>>
//!\endcond
class trace_matrix_simd_adapter_iterator :
    public two_dimensional_matrix_iterator_base<trace_matrix_simd_adapter_iterator<matrix_iter_t>,
                                                matrix_major_order::column>
{
private:
    //!\brief The base class.
    using base_t = two_dimensional_matrix_iterator_base<trace_matrix_simd_adapter_iterator<matrix_iter_t>,
                                                        matrix_major_order::column>;
    //!\brief Befriend the base class.
    friend base_t;

    //!\brief Befriend the corresponding const iterator.
    template <two_dimensional_matrix_iterator other_matrix_iter_t>
        requires simd::simd_concept<std::iter_value_t<other_matrix_iter_t>>
    friend class trace_matrix_simd_adapter_iterator;

    matrix_iter_t host_iter{}; //!< The wrapped iterator.
    size_t simd_lane{}; //!< The simd vector lane to get the seqan3::detail::trace_directions from.

public:
    /*!\name Associated types
    * \{
    */
    //!\brief Value type of this iterator.
    using value_type = trace_directions;
    //!\brief Reference is the same as the `value_type`.
    using reference = trace_directions;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Type for distances between iterators.
    using difference_type = std::iter_difference_t<matrix_iter_t>;
    //!\brief The iterator tag.
    using iterator_category = std::random_access_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
    * \{
    */
    constexpr trace_matrix_simd_adapter_iterator() = default; //!< Defaulted.
    constexpr trace_matrix_simd_adapter_iterator(trace_matrix_simd_adapter_iterator const &)
        = default; //!< Defaulted.
    constexpr trace_matrix_simd_adapter_iterator(trace_matrix_simd_adapter_iterator &&)
        = default; //!< Defaulted.
    constexpr trace_matrix_simd_adapter_iterator & operator=(trace_matrix_simd_adapter_iterator const &)
        = default; //!< Defaulted.
    constexpr trace_matrix_simd_adapter_iterator & operator=(trace_matrix_simd_adapter_iterator &&)
        = default; //!< Defaulted.
    ~trace_matrix_simd_adapter_iterator() = default; //!< Defaulted.

    /*!\brief Construction from the underlying matrix iterator and the corresponding simd lane to represent.
     *
     * \param[in] iter The underlying iterator over the actual matrix.
     * \param[in] simd_lane The simd lane to read from.
     */
    constexpr trace_matrix_simd_adapter_iterator(matrix_iter_t iter, size_t simd_lane) :
        host_iter{std::move(iter)},
        simd_lane{std::move(simd_lane)}
    {
        assert(this->simd_lane < simd::simd_traits<std::iter_value_t<matrix_iter_t>>::length);
    }

    //!\brief Construction of cons-iterator from non-const-iterator.
    template <typename other_matrix_iter_t>
        requires (!std::same_as<other_matrix_iter_t, matrix_iter_t> &&
                  std::constructible_from<matrix_iter_t, other_matrix_iter_t>)
    constexpr trace_matrix_simd_adapter_iterator(trace_matrix_simd_adapter_iterator<other_matrix_iter_t> other) noexcept
        : host_iter{std::move(other.host_iter)},
          simd_lane{std::move(other.simd_lane)}
    {}
    //!\}

    //!\brief Dereferences the wrapped iterator and converts the value at `simd_lane` to a
    //!\      seqan3::detail::trace_directions.
    constexpr reference operator*() const noexcept
    {
        assert(is_valid_trace_directions((*host_iter)[simd_lane]));

        using enum_value_t = std::underlying_type_t<trace_directions>;
        return trace_directions{enum_value_t((*host_iter)[simd_lane])};
    }

    // Import operator+= of base class.
    using base_t::operator+=;

    //!\brief Advances the iterator by the given `offset`.
    constexpr trace_matrix_simd_adapter_iterator & operator+=(matrix_offset const & offset) noexcept
    {
        host_iter += offset;
        return *this;
    }

    //!\copydoc seqan3::detail::two_dimensional_matrix_iterator::coordinate()
    constexpr matrix_coordinate coordinate() const noexcept
    {
        return host_iter.coordinate();
    }

private:
    /*!\brief Helper to check if the stored simd value is a valid seqan3::trace_directions.
     * \param[in] direction The direction to check.
     *
     * \returns `true` if it is a seqan3::trace_directions, otherwise `false`.
     */
    constexpr bool is_valid_trace_directions(int32_t const direction) const noexcept
    {
        return (0 <= direction) && (direction < 32);
    }
};
}  // namespace seqan3::detail
