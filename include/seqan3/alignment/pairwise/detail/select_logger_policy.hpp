// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::select_logger_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger_dev_null.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger.hpp>
#include <seqan3/core/detail/empty_type.hpp>

namespace seqan3::detail
{

/*!\brief Selects the algorithm logger policy from the algorithm type traits.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam algorithm_traits_t The type of the configured algorithm traits.
 *
 * \details
 *
 * Configures the algorithm logger policy type based on the configured alignment algorithm traits.
 * The selector differentiates between a standard logger that stores each column of the alignment matrix or
 * a dev_null logger, which offers the interfaces but does not log any information, as if by sending output to
 * /dev/null. The latter one is chosen if the alignment is not configured with debug information.
 *
 * If the debug information is enabled but the trace is not computed, then only the debug score matrix is configured
 * and the debug trace matrix defaults to seqan3::detail::empty_type.
 */
template <typename algorithm_traits_t>
struct select_logger_policy
{
private:
    //!\brief The configured score type for logging.
    using debug_score_t = std::optional<typename algorithm_traits_t::score_type>;
    //!\brief The configured trace type for logging.
    using debug_trace_t = std::optional<typename algorithm_traits_t::trace_type>;
    //!\brief The debug score matrix type.
    using debug_score_matrix_t = two_dimensional_matrix<debug_score_t,
                                                        std::allocator<debug_score_t>,
                                                        matrix_major_order::column>;

    //!\brief The debug trace matrix type.
    using debug_trace_matrix_t = std::conditional_t<algorithm_traits_t::compute_sequence_alignment,
                                                    two_dimensional_matrix<debug_trace_t,
                                                                           std::allocator<debug_trace_t>,
                                                                           matrix_major_order::column>,
                                                    empty_type>;
public:
    //!\brief The configured logger policy.
    using type = std::conditional_t<algorithm_traits_t::is_debug,
                                    policy_alignment_algorithm_logger<debug_score_matrix_t, debug_trace_matrix_t>,
                                    policy_alignment_algorithm_logger_dev_null>;
};

/*!\brief Helper type alias for selecting the type of the algorithm logger policy.
 * \ingroup pairwise_alignment
 * \relates seqan3::detail::select_logger_policy
 */
template <typename algorithm_traits_t>
using select_logger_policy_t = typename select_logger_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
