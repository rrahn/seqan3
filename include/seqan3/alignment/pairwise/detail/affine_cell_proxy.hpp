// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_cell_proxy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/core/tuple_utility.hpp>

namespace seqan3::detail
{

template <tuple_like tuple_t>
class affine_cell_proxy : public tuple_t
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    affine_cell_proxy() = default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy &&) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy &&) = default; //!< Defaulted.
    ~affine_cell_proxy() = default; //!< Defaulted.

    using tuple_t::tuple_t;
    using tuple_t::operator=;
    //!\}

    decltype(auto) optimal_score() & { return std::get<0>(*this); }
    decltype(auto) optimal_score() const & { return std::get<0>(*this); }
    decltype(auto) optimal_score() && { return std::get<0>(std::move(*this)); }
    decltype(auto) optimal_score() const && { return std::get<0>(std::move(*this)); }

    decltype(auto) horizontal_score() & { return std::get<1>(*this); }
    decltype(auto) horizontal_score() const & { return std::get<1>(*this); }
    decltype(auto) horizontal_score() && { return std::get<1>(std::move(*this)); }
    decltype(auto) horizontal_score() const && { return std::get<1>(std::move(*this)); }

    decltype(auto) vertical_score() & { return std::get<2>(*this); }
    decltype(auto) vertical_score() const & { return std::get<2>(*this); }
    decltype(auto) vertical_score() const && { return std::get<2>(std::move(*this)); }
};
} // namespace seqan3::detail

namespace std
{

template <typename tuple_t>
struct tuple_size<seqan3::detail::affine_cell_proxy<tuple_t>> : tuple_size<tuple_t>
{};

template <size_t index, typename tuple_t>
struct tuple_element<index, seqan3::detail::affine_cell_proxy<tuple_t>> : tuple_element<index, tuple_t>
{};
} // namespace std

