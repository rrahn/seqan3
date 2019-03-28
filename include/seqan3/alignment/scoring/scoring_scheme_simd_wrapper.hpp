// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::scoring_scheme_simd_wrapper.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace seqan3::detail
{

template <Simd simd_t>
class scoring_scheme_simd_wrapper
{
private:

    //!\brief The alternative simd scoring schemes.
    using scheme_variant_t = std::variant<simd_scoring_scheme_simple<simd_t>>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr scoring_scheme_simd_wrapper()                                                = default; //!< Defaulted.
    constexpr scoring_scheme_simd_wrapper(scoring_scheme_simd_wrapper const &)             = default; //!< Defaulted.
    constexpr scoring_scheme_simd_wrapper(scoring_scheme_simd_wrapper &&)                  = default; //!< Defaulted.
    constexpr scoring_scheme_simd_wrapper & operator=(scoring_scheme_simd_wrapper const &) = default; //!< Defaulted.
    constexpr scoring_scheme_simd_wrapper & operator=(scoring_scheme_simd_wrapper &&)      = default; //!< Defaulted.
    ~scoring_scheme_simd_wrapper()                                                         = default; //!< Defaulted.

    template <typename scoring_scheme_t>
    constexpr scoring_scheme_simd_wrapper(scoring_scheme_t && scheme)
    {
        // Parse the scheme and decide which scheme to use.

        schemes = simd_scoring_scheme_simple{std::forward<scoring_scheme_t>(scheme)};
    }
    //!\}

    /*!\brief Computes the score of two sequence entries.
     * \param lhs The query simd vector.
     * \param rhs The database simd vector.
     */
    template <typename alignment_t>
    constexpr simd_t score(simd_t const & lhs, simd_t const & rhs)
    {
        return std::visit([&](auto & scheme) { return scheme.score<alignment_t>(lhs, rhs); }, schemes);
    }

private:

    //!\brief Holds the alternative simd scoring schemes.
    scheme_variant_t schemes;
};

} // namespace seqan3::detail
