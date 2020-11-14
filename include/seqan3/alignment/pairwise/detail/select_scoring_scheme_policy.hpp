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
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?

template <typename algorithm_traits_t, typename alphabet_t, bool is_vectorised>
struct select_scoring_scheme_policy
{
public:
    using type = policy_scoring_scheme<typename algorithm_traits_t::scoring_scheme_type>;
};

template <typename algorithm_traits_t, typename alphabet_t>
struct select_scoring_scheme_policy<algorithm_traits_t, alphabet_t, true>
{
private:
    using alignment_method_t = std::conditional_t<algorithm_traits_t::is_local,
                                                  seqan3::align_cfg::method_local,
                                                  seqan3::align_cfg::method_global>;

    static constexpr bool is_aminoacid_scheme =
        is_type_specialisation_of_v<typename algorithm_traits_t::scoring_scheme_type,
                                    aminoacid_scoring_scheme>;

    using simple_simd_scheme_t = simd_match_mismatch_scoring_scheme<typename algorithm_traits_t::score_type,
                                                                    alphabet_t,
                                                                    alignment_method_t>;

    using matrix_simd_scheme_t = simd_matrix_scoring_scheme<typename algorithm_traits_t::score_type,
                                                            alphabet_t,
                                                            alignment_method_t>;

    using alignment_scoring_scheme_t = std::conditional_t<is_aminoacid_scheme,
                                                          matrix_simd_scheme_t,
                                                          simple_simd_scheme_t>;
public:
    using type = policy_scoring_scheme<alignment_scoring_scheme_t>;
};

template <typename algorithm_traits_t, typename alphabet_t>
using select_scoring_scheme_policy_t =
    typename select_scoring_scheme_policy<algorithm_traits_t,
                                          alphabet_t,
                                          algorithm_traits_t::is_vectorised>::type;

}  // namespace seqan3::detail
