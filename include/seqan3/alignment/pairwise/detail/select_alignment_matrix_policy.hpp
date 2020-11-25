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

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_affine_cell.hpp>
#include <seqan3/alignment/matrix/detail/combined_score_and_trace_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_matrix.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?

template <typename algorithm_traits_t>
struct select_alignment_matrix_policy
{
private:

    using matrix_element_t = alignment_matrix_element_affine_cell<typename algorithm_traits_t::score_type>;
    using score_matrix_t = score_matrix_single_column<matrix_element_t>;
    using trace_matrix_t = trace_matrix_full<typename algorithm_traits_t::trace_type>;

    using alignment_matrix_t = std::conditional_t<algorithm_traits_t::requires_trace_information,
                                                  combined_score_and_trace_matrix<score_matrix_t, trace_matrix_t>,
                                                  score_matrix_t>;

    using alignment_matrix_traits_t = alignment_matrix_traits<typename algorithm_traits_t::score_type,
                                                              typename algorithm_traits_t::matrix_index_type,
                                                              algorithm_traits_t::requires_trace_information>;

public:
    //!\brief The configured recursion policy.
    using type = policy_alignment_matrix<alignment_matrix_traits_t, alignment_matrix_t>;
};

template <typename algorithm_traits_t>
using select_alignment_matrix_policy_t = typename select_alignment_matrix_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
