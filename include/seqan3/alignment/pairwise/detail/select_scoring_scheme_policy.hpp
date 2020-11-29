// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_scoring_scheme_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/detail/simd_matrix_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/pairwise/detail/policy_scoring_scheme.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{

/*!\brief Selects the scoring scheme policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 * \tparam alphabet_t The type of the alphabet to generate the simd alignment scoring scheme for; must model
 *                    seqan3::semialphabet.
 */
template <typename algorithm_traits_t, semialphabet alphabet_t>
struct select_scoring_scheme_policy
{
private:
    //!\brief A flag indicating whether an amino acid scoring scheme was selected.
    static constexpr bool is_aminoacid_scheme =
        is_type_specialisation_of_v<typename algorithm_traits_t::scoring_scheme_type, aminoacid_scoring_scheme>;
    //!\brief The configured score type.
    using score_t = typename algorithm_traits_t::score_type;
    //!\brief The configured alignment method type.
    using method_t = std::conditional_t<algorithm_traits_t::is_local,
                                        align_cfg::method_local,
                                        align_cfg::method_global>;

    //!\brief The selected simd scoring scheme.
    template <typename simd_score_t>
    using simd_scoring_scheme_t =
        std::conditional_t<is_aminoacid_scheme,
                           simd_matrix_scoring_scheme<simd_score_t, alphabet_t, method_t>,
                           simd_match_mismatch_scoring_scheme<simd_score_t, alphabet_t, method_t>>;

    //!\brief The selected scoring scheme type.
    using scoring_scheme_t = lazy_conditional_t<algorithm_traits_t::is_vectorised,
                                                lazy<simd_scoring_scheme_t, score_t>,
                                                typename algorithm_traits_t::scoring_scheme_type>;
public:
    //!\brief The selected scoing scheme policy.
    using type = policy_scoring_scheme<scoring_scheme_t>;
};

/*!\brief Helper type alias for selecting the type of the scoring scheme policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_scoring_scheme_policy
 */
template <typename algorithm_traits_t, semialphabet alphabet_t>
using select_scoring_scheme_policy_t = typename select_scoring_scheme_policy<algorithm_traits_t, alphabet_t>::type;

}  // namespace seqan3::detail
