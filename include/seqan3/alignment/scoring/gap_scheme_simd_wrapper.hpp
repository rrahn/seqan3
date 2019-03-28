// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::simd_gap_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

namespace seqan3::detail
{

// ------------------------------------------------------------------
// seqan3::simd_gap_scheme
// ------------------------------------------------------------------

/*!\brief A scheme for representing and computing scores against gap characters.
 * \tparam score_type Type of the score values saved internally.
 * \ingroup scoring
 */
template <Simd simd_t>
class gap_scheme_simd_wrapper
{
public:

    //!\brief The score type.
    using score_type = simd_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Initialises the gap scheme with -1 and 0 for gap and gap open score respectively.
    constexpr gap_scheme_simd_wrapper() : gap_score{simd::fill<simd_t>(-1)}, gap_open_score{simd::fill<simd_t>(0)}
    {}
    constexpr gap_scheme_simd_wrapper(gap_scheme_simd_wrapper const &)             noexcept = default; //!< Defaulted
    constexpr gap_scheme_simd_wrapper(gap_scheme_simd_wrapper &&)                  noexcept = default; //!< Defaulted
    constexpr gap_scheme_simd_wrapper & operator=(gap_scheme_simd_wrapper const &) noexcept = default; //!< Defaulted
    constexpr gap_scheme_simd_wrapper & operator=(gap_scheme_simd_wrapper &&)      noexcept = default; //!< Defaulted
    ~gap_scheme_simd_wrapper()                                                     noexcept = default; //!< Defaulted

    /*!\brief Constructs from the seqan3::gap_scheme.
     */
    template <Arithmetic score_type>
    constexpr gap_scheme_simd_wrapper(gap_scheme<score_type> const & gap_scheme)
    {
        init(gap_scheme);
    }

    template <Arithmetic score_type>
    constexpr gap_scheme_simd_wrapper(gap_scheme<score_type> && gap_scheme)
    {
        init(std::move(gap_scheme));
    }
    //!\}

    /*!\name Setter
     * \{
     */
    template <Arithmetic score_type>
    constexpr void set_gap_scheme(gap_scheme<score_type> const & gap_scheme)
    {
        init(gap_scheme);
    }

    template <Arithmetic score_type>
    constexpr void set_gap_scheme(gap_scheme<score_type> && gap_scheme)
    {
        init(std::move(gap_scheme));
    }
    //!\}

    /*!\name Gap accessors
     * \{
     */

    constexpr simd_t & get_gap_score() noexcept
    {
        return gap_score;
    }

    constexpr simd_t const & get_gap_score() const noexcept
    {
        return gap_score;
    }

    constexpr simd_t & get_gap_open_score() noexcept
    {
        return gap_open_score;
    }

    constexpr simd_t const & get_gap_open_score() const noexcept
    {
        return gap_open_score;
    }
    //!\}

private:

    template <typename gap_scheme_t>
    void init(gap_scheme_t && gap_scheme)
    {
        using score_type  = typename std::remove_reference_t<gap_scheme_t>::score_type;
        using scalar_type = typename simd_traits<simd_t>::scalar_type;

        // Check if scalar type and score type are compatible.
        if constexpr (FloatingPoint<score_type>)
            static_assert(FloatingPoint<scalar_type>, "The simd type is not based on a floating point type.");
        else
            static_assert(std::Integral<scalar_type>, "The simd type is not based on an integral type.");

        if constexpr (sizeof(scalar_type) < sizeof(score_type))
        {
            scalar_type min_value = std::numeric_limits<scalar_type>::lowest();
            if (gap_scheme.get_gap_score() < static_cast<score_type>(min_value))
                throw std::invalid_argument{"The selected gap score overflows with the current simd type."};
            else if (gap_scheme.get_gap_open_score() < static_cast<score_type>(min_value))
                throw std::invalid_argument{"The selected gap open score overflows with the current simd type."};
        }

        gap_score      = simd::fill<simd_t>(gap_scheme.get_gap_score());
        gap_open_score = simd::fill<simd_t>(gap_scheme.get_gap_open_score());
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::CerealArchive.
     * \param  archive   The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <CerealArchive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(gap_score);
        archive(gap_open_score);
    }
    //!\endcond

    //!\brief The gap score.
    simd_t gap_score;
    //!\brief The gap open score.
    simd_t gap_open_score;
};

/*!\name Type deduction guides
 * \relates seqan3::gap_scheme_simd_wrapper
 * \{
 */

/*!\brief This guide deduces the simd type from the score type of the passed seqan3::gap_scheme.
 * To use a custom simd type, specify the template argument manually.
 */
template <typename score_type>
gap_scheme_simd_wrapper(gap_scheme<score_type>) -> gap_scheme_simd_wrapper<typename simd_type<score_type>::type>;

//!\}

} // namespace seqan3::detail
