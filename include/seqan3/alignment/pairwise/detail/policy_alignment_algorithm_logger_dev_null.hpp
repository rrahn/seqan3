// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_algorithm_logger_dev_null.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Implements a /dev/null-like logger by ignoring the data.
 * \ingroup pairwise_alignment
 *
 * \details
 *
 * This logger ignores all input that is send to it. It does not create any data members and calling the respective
 * interfaces will be optimised away by the compiler.
 */
class policy_alignment_algorithm_logger_dev_null
{
protected:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_algorithm_logger_dev_null() = default; //!< Defaulted.
    policy_alignment_algorithm_logger_dev_null(policy_alignment_algorithm_logger_dev_null const &)
        = default; //!< Defaulted.
    policy_alignment_algorithm_logger_dev_null(policy_alignment_algorithm_logger_dev_null &&)
        = default; //!< Defaulted.
    policy_alignment_algorithm_logger_dev_null & operator=(policy_alignment_algorithm_logger_dev_null const &)
        = default; //!< Defaulted.
    policy_alignment_algorithm_logger_dev_null & operator=(policy_alignment_algorithm_logger_dev_null &&)
        = default; //!< Defaulted.
    ~policy_alignment_algorithm_logger_dev_null() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_alignment_algorithm_logger::policy_alignment_algorithm_logger(alignment_configuration_t const &)
    template <typename alignment_configuration_t>
    policy_alignment_algorithm_logger_dev_null(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {}
    //!\}

    //!\copydoc seqan3::detail::policy_alignment_algorithm_logger::initialise_debug_matrices
    void initialise_debug_matrices(size_t const SEQAN3_DOXYGEN_ONLY(sequence1_size),
                                   size_t const SEQAN3_DOXYGEN_ONLY(sequence2_size))
    {}

    //!\copydoc seqan3::detail::policy_alignment_algorithm_logger::log_alignment_matrix_column
    template <typename coordinate_column_t, typename alignment_column_t>
    void log_alignment_matrix_column(coordinate_column_t && SEQAN3_DOXYGEN_ONLY(coordinate_column),
                                     alignment_column_t && SEQAN3_DOXYGEN_ONLY(alignment_column))
    {}
};

}  // namespace seqan3::detail
