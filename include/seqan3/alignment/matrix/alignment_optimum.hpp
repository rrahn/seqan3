// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_optimum.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief Stores the current optimum of the alignment algorithms.
 * \ingroup alignment_matrix
 * \tparam score_t The type of the tracked alignment score; must model seqan3::Arithmetic or seqan3::Simd.
 *
 * \details
 *
 * This is an aggregate type, so the score needs to be passed before the seqan3::alignment_coordinate during
 * construction.
 */
template <typename score_t>
//!\cond
    requires Arithmetic<score_t> || Simd<score_t>
//!\endcond
struct alignment_optimum
{
    //!\brief The optimal score.
    score_t score;
    //!\brief The corresponding coordinate within the alignment matrix.
    alignment_coordinate<score_t> coordinate;

    constexpr friend alignment_optimum<score_t> max(alignment_optimum const & lhs, alignment_optimum const & rhs) noexcept
    {
        if constexpr (Simd<score_t>)
        {
            alignment_optimum<score_t> tmp{};
            auto mask = lhs.score < rhs.score;
            tmp.score = (mask) ? rhs.score : lhs.score;
            tmp.coordinate.first = (mask) ? rhs.coordinate.first : lhs.coordinate.first;
            tmp.coordinate.second = (mask) ? rhs.coordinate.second : lhs.coordinate.second;

            return tmp;
        }
        else
        {
            return (lhs.score < rhs.score) ? rhs : lhs;
        }
    }
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Default constructed objects deduce to `int32_t`.
alignment_optimum() -> alignment_optimum<int32_t>;

//!\brief Deduce the score type.
template <typename score_t>
//!\cond
    requires Arithmetic<score_t> || Simd<score_t>
//!\endcond
alignment_optimum(score_t const &, alignment_coordinate<score_t> const &) -> alignment_optimum<score_t>;
//!\}

/*!\brief A less than comparator for two seqan3::detail::alignment_optimum objects.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This function object is used in std::max functions to compare two seqan3::detail::alignment_optimum objects.
 */
struct alignment_optimum_compare_less
{

    /*!\brief Function call operator that implements less than comparison.
     * \tparam lhs_t The type of the left-hand side operand. Must be a type specialisation of
     *               seqan3::detail::alignment_optimum.
     * \tparam rhs_t The type of the right-hand side operand. Must be a type specialisation of
     *               seqan3::detail::alignment_optimum.
     *
     * \param[in] lhs The left-hand side operand.
     * \param[in] rhs The right-hand side operand.
     *
     * \returns bool `true` if `lhs.score < rhs.score`, otherwise `false`.
     */
    template <typename lhs_t, typename rhs_t>
    //!\cond
        requires (is_type_specialisation_of_v<lhs_t, alignment_optimum> &&
                  is_type_specialisation_of_v<rhs_t, alignment_optimum>)
    //!\endcond
    constexpr bool operator()(lhs_t const & lhs, rhs_t const & rhs) const
    {
        return lhs.score < rhs.score;
    }

};

} // namespace seqan3::detail
