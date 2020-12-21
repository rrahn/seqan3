// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::static_thread_pool_operation_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::execution
{

class static_thread_pool_operation_base
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr static_thread_pool_operation_base() = default; //!< Defaulted.
    constexpr static_thread_pool_operation_base(static_thread_pool_operation_base const &) = delete; //!< Deleted.
    constexpr static_thread_pool_operation_base(static_thread_pool_operation_base &&) = delete; //!< Deleted.
    constexpr static_thread_pool_operation_base & operator=(static_thread_pool_operation_base const &) = delete; //!< Deleted.
    constexpr static_thread_pool_operation_base & operator=(static_thread_pool_operation_base &&) = delete; //!< Deleted.
    ~static_thread_pool_operation_base() = default; //!< Defaulted.
    //!\}

    //!\brief Invokes the concrete operation.
    virtual void set_value() noexcept = 0;
    virtual void set_error() noexcept = 0;
    virtual void set_done() noexcept = 0;
};

}  // namespace seqan3::execution
