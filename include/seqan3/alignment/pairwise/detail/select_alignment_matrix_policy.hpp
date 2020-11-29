// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_alignment_matrix_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/combined_score_and_trace_matrix.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_matrix.hpp>

namespace seqan3::detail
{

/*!\brief Selects the alignment matrix policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 */
template <typename algorithm_traits_t>
struct select_alignment_matrix_policy
{
private:
    //!\brief The type of the score matrix.
    using score_matrix_t = score_matrix_single_column<typename algorithm_traits_t::score_type>;
    //!\brief The type of the trace matrix.
    using trace_matrix_t = trace_matrix_full<typename algorithm_traits_t::trace_type>;
    //!\brief The type of the alignment matrix.
    using alignment_matrix_t = std::conditional_t<algorithm_traits_t::requires_trace_information,
                                                  combined_score_and_trace_matrix<score_matrix_t, trace_matrix_t>,
                                                  score_matrix_t>;
    //!\brief The type of the coordinate matrix.
    using coordinate_matrix_t = coordinate_matrix<typename algorithm_traits_t::matrix_index_type>;
public:
    //!\brief The type of the selected alignment matrix policy.
    using type = policy_alignment_matrix<alignment_matrix_t, coordinate_matrix_t>;
};

/*!\brief Helper type alias for selecting the type of the alignment matrix policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_alignment_matrix_policy
 */
template <typename algorithm_traits_t>
using select_alignment_matrix_policy_t = typename select_alignment_matrix_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
