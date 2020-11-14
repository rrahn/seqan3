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
#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger.hpp>
#include <seqan3/core/detail/empty_type.hpp>

namespace seqan3::detail
{
// TODO: Can we also work with forward declarations here?

class empty_policy
{
public:
    empty_policy() = default;
    template <typename configuration_t>
    empty_policy(configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {}

    void initialise_debug_matrices(size_t const SEQAN3_DOXYGEN_ONLY(sequence1_size),
                                   size_t const SEQAN3_DOXYGEN_ONLY(sequence2_size))
    {}

    template <typename coordinate_column_t, typename alignment_column_t>
    void log_alignment_matrix_column(coordinate_column_t && SEQAN3_DOXYGEN_ONLY(coordinate_column),
                                     alignment_column_t && SEQAN3_DOXYGEN_ONLY(alignment_column))
    {}
};

template <typename algorithm_traits_t>
struct select_logger_policy
{
private:

    using debug_score_t = std::optional<typename algorithm_traits_t::score_type>;
    using debug_trace_t = std::optional<typename algorithm_traits_t::trace_type>;
    using debug_score_matrix_t = two_dimensional_matrix<debug_score_t,
                                                        std::allocator<debug_score_t>,
                                                        matrix_major_order::column>;

    using debug_trace_matrix_t = std::conditional_t<algorithm_traits_t::compute_sequence_alignment,
                                                    two_dimensional_matrix<debug_trace_t,
                                                                           std::allocator<debug_trace_t>,
                                                                           matrix_major_order::column>,
                                                    empty_type>;
public:
    //!\brief The configured recursion policy.
    using type = std::conditional_t<algorithm_traits_t::is_debug,
                                    policy_alignment_algorithm_logger<debug_score_matrix_t, debug_trace_matrix_t>,
                                    empty_policy>;
};

template <typename algorithm_traits_t>
using select_logger_policy_t = typename select_logger_policy<algorithm_traits_t>::type;

}  // namespace seqan3::detail
