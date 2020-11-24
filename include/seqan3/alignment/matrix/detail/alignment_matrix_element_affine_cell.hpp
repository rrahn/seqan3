// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_element_affine_cell.hpp.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <tuple>
#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

template <typename tuple_t>
class affine_cell : public tuple_t
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    affine_cell() = default; //!< Defaulted.
    affine_cell(affine_cell const &) = default; //!< Defaulted.
    affine_cell(affine_cell &&) = default; //!< Defaulted.
    affine_cell & operator=(affine_cell const &) = default; //!< Defaulted.
    affine_cell & operator=(affine_cell &&) = default; //!< Defaulted.
    ~affine_cell() = default; //!< Defaulted.

    // Inherit the base class's constructor to enable element-wise initialisation (direct and converting constructor).
    using tuple_t::tuple_t;

    //!\brief Converting copy-constructor.
    template <typename other_tuple_t>
    //!\cond
        requires std::constructible_from<tuple_t, other_tuple_t &&>
    //!\endcond
    explicit affine_cell(affine_cell<other_tuple_t> other) :
        tuple_t{static_cast<other_tuple_t &&>(other)}
    {}

    //!\brief Converting copy-assignment.
    template <typename other_tuple_t>
    //!\cond
        requires std::assignable_from<tuple_t &, other_tuple_t &&>
    //!\endcond
    affine_cell & operator=(affine_cell<other_tuple_t> other)
    {
        as_base() = static_cast<other_tuple_t &&>(other);
        return *this;
    }
    //!\}

    /*!\name Score value accessor
     * \brief Specific accessor function to get the respective score value from an affine matrix cell.
     * \{
     */
    //!\brief Access the best score of the wrapped score matrix cell.
    decltype(auto) best_score() & noexcept { return get_score_impl<0>(*this); }
    //!\overload
    decltype(auto) best_score() const & noexcept { return get_score_impl<0>(*this); }
    //!\overload
    decltype(auto) best_score() && noexcept { return get_score_impl<0>(std::move(*this)); }

    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_score() & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() const & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() && noexcept { return get_score_impl<1>(std::move(*this)); }

    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_score() & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() const & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() && noexcept { return get_score_impl<2>(std::move(*this)); }
    //!\}

private:
    /*!\brief Implements the get interface for the various calls to receive the score value.
     * \tparam index The index of the tuple element to get; must be smaller than 3.
     * \tparam this_t The perfectly forwarded type of `*this`.
     *
     * \param[in] me The instance of `*this`.
     *
     * \returns The score value from the given tuple index.
     */
    template <size_t index, typename this_t>
    //!\cond
        requires (index < 3)
    //!\endcond
    static constexpr decltype(auto) get_score_impl(this_t && me) noexcept
    {
        using std::get;

        return get<index>(std::forward<this_t>(me));
    }

    //!\brief Casts `this` to the base class type.
    tuple_t & as_base() & noexcept
    {
        return static_cast<tuple_t &>(*this);
    }
};

template <std::semiregular score_t>
struct alignment_matrix_element_affine_cell
{
private:
    score_t vertical_score{};
public:

    using value_type = score_t;
    using storage_value_type = std::pair<value_type, value_type>;
    using element_type = affine_cell<std::tuple<value_type, value_type, value_type>>;
    using element_reference = affine_cell<std::tuple<value_type &, value_type &, value_type &>>;

    static constexpr bool returns_element_proxy = true;

    constexpr element_reference make_element(storage_value_type & storage_value)
    {
        return {storage_value.first, storage_value.second, vertical_score};
    }

    constexpr storage_value_type initialise(value_type const value) noexcept
    {
        vertical_score = value;
        return {value, value};
    }
};
}  // namespace seqan3::detail
