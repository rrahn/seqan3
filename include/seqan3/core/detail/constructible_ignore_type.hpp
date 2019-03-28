// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::constructible_ignore_type.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace seqan3::detail
{

class constructible_ignore_type
{
public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr constructible_ignore_type()                                              = default; //!< Defaulted.
    constexpr constructible_ignore_type(constructible_ignore_type const &)             = default; //!< Defaulted.
    constexpr constructible_ignore_type(constructible_ignore_type &&)                  = default; //!< Defaulted.
    constexpr constructible_ignore_type & operator=(constructible_ignore_type const &) = default; //!< Defaulted.
    constexpr constructible_ignore_type & operator=(constructible_ignore_type &&)      = default; //!< Defaulted.
    ~constructible_ignore_type()                                                       = default; //!< Defaulted.

    /*!\brief Ignores constructor argument.
     * \tparam t The type of the argument to ignore.
     */
    template <typename t>
    constexpr constructible_ignore_type(t const &) noexcept
    {}

    /*!\brief Ignores the assigned argument.
     * \tparam t The type of the assigned argument.
     */
    template <typename t>
    constexpr constructible_ignore_type & operator=(t const &) noexcept
    {
        return *this;
    }
    //!\}
};

} // namespace seqan3::detail
