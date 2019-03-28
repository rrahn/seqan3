// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::simd_scoring_scheme_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/simd/all.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief A simple scoring scheme for vectorised alignments.
 * \ingroup scoring
 * \tparam simd_t       The type of the SIMD vector.
 * \tparam alignment_t  The type of the alignment to compute.
 *
 * \details
 *
 * Provides a scoring scheme that compares packed elements in two SIMD vectors and returns a new SIMD vector packed
 * with match score or mismatch score depending on the result of the comparison.
 */
template <Simd simd_t, Alphabet alphabet_t>
class simd_scoring_scheme_simple
{
public:

    // We require alphabets with at least 2 values.
    static_assert(alphabet_size_v<alphabet_t> > 1, "Only alphabets with at least 2 values are allowed.");

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_scoring_scheme_simple() = default;                                                  //!< Default.
    constexpr simd_scoring_scheme_simple(simd_scoring_scheme_simple const &) = default;                //!< Default.
    constexpr simd_scoring_scheme_simple(simd_scoring_scheme_simple &&) = default;                     //!< Default.
    constexpr simd_scoring_scheme_simple & operator=(simd_scoring_scheme_simple const &) = default;    //!< Default.
    constexpr simd_scoring_scheme_simple & operator=(simd_scoring_scheme_simple &&) = default;         //!< Default.
    ~simd_scoring_scheme_simple() = default;                                                           //!< Default.

    /*!\brief Creates a new scoring scheme from the passed match and mismatch score.
     * \tparam  score_type The type of the match and mismatch score. Must be convertible to the scalar type of `simd_t`.
     * \param[in] match    The match score.
     * \param[in] mismatch The mismatch score.
     */
    template <ScoringScheme<alphabet_t> scoring_scheme_t>
    constexpr simd_scoring_scheme_simple(scoring_scheme_t && scheme)
    {
        using score_t = typename scoring_scheme_t::score_type;
        using scalar_t = typename simd_traits<simd_t>::scalar_type;

        if constexpr (FloatingPoint<score_t>)
            static_assert(FloatingPoint<scalar_t>, "The simd type is not based on floating point types.");
        else
            static_assert(std::Integral<scalar_t>, "The simd type is not based on integral types.");

        if constexpr (sizeof(scalar_t) < sizeof(score_t))
        {
            if (scheme.score(assign_rank_to(0, alphabet_t{}), assign_rank_to(0, alphabet_t{})) >
                static_cast<score_t>(std::numeric_limits<scalar_t>::max()) ||
                scheme.score(assign_rank_to(0, alphabet_t{}), assign_rank_to(1, alphabet_t{})) <
                static_cast<score_t>(std::numeric_limits<scalar_t>::lowest()))
            {
                throw std::invalid_argument{"The selected score scheme score overflows with the current simd type."};
            }
        }

        match_score = simd::fill<simd_t>(scheme.score(assign_rank_to(0, alphabet_t{}),
                                                      assign_rank_to(0, alphabet_t{})));
        mismatch_score = simd::fill<simd_t>(scheme.score(assign_rank_to(0, alphabet_t{}),
                                                         assign_rank_to(1, alphabet_t{})));
    }
    //!\}

    /*!\brief Computes the score of two SIMD vectors.
     * \tparam alignment_t The alignment type used for the comparison.
     * \param[in] lhs The left operand to compare.
     * \param[in] rhs The right operand to compare.
     * \returns
     *
     * \details
     *
     * This function compares packed elements in both vectors and returns a new vector filled with match and mismatch
     * scores depending on the result of the comparison.
     * For sequences with different size, a special padding sequence is used. The sequences are padded in such a way
     * that comparison with the padding characters either produces a mismatch or a match.
     * In case of the local alignment the vertical sequences are padded with a value that is distinct to the any letter
     * of the corresponding alphabet and to the padding character of the horizontal sequences.
     * Thus XOR-ing vectors from the vertical and horizontal sequences will only result in 0 if the characters are
     * identical and none of the characters (horizontal or vertical) are padding characters.
     * For the global alignment we want all padding characters to always produce a match in order to propagate the
     * correct score to the end of the largest matrix computation. In this case we artificially set the signed bit of
     * the packed vector elements and test whether the XOR result is less or equal to 0.
     * This is only true if both characters are identical or if one of them is a padding character. For the local
     * alignment we want to explicitly produce mismatches with padded characters.
     *
     * ### Exception
     *
     * no-throw guarantee.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    template <typename alignment_t>
    constexpr simd_t score(simd_t const & lhs, simd_t const & rhs) const noexcept
    {
        static_assert(std::Same<alignment_t, detail::local_alignment_type> ||
                      std::Same<alignment_t, detail::global_alignment_type>,
                      "This alignment type is not supported.");

        typename simd_traits<simd_t>::mask_type mask;
        // For global and local alignment we have slightly different formulas.
        // In global alignment we want padding characters always to match
        if constexpr (std::Same<alignment_t, detail::global_alignment_type>)
            mask = (lhs ^ rhs) <= 0;
        else // local alignment type; padded characters should always mismatch.
            mask = (lhs ^ rhs) == 0;

        return mask ? match_score : mismatch_score;

        // global alignment:
        // lhs ^ rhs &

        // out of the alignment produce match.
        // a == f -> match

        // local alignment: last two bits set for padded chars.
        // global alignment: all elements added with last bit.

        // regular char -> comparison
        // padded char -> bit is set in last element

        // in global alignment same mask
        // 1100 ^ 1101 = 0000  -> 1 match    ->  0000 & 0001 = 0
        // 1100 ^ 1001 = 0100  -> 0 mismatch ->  0100 & 0001 = 0
        // 1100 ^ 0001 = 1101  -> 0 mismatch ->  1101 & 0001 = needs to be a match.
        // 0000 ^ 0001 = 0001  -> 1 match
        // 0001 ^ 0001 = 0000  -> 1 match

        // 0001 ^ 0010 = 0011  -> mismatch
        // 0001 ^ 0001 = 0000  ->match (not in local alignment)
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
        archive(match_score);
        archive(mismatch_score);
    }
    //!\endcond

private:

    //!\brief The score for a match.
    simd_t match_score;
    //!\brief The score for a mismatch.
    simd_t mismatch_score;
};

} // namespace seqan3::detail
