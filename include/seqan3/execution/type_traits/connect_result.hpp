// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides type traits for the execution.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <seqan3/execution/cpo/connect.hpp>

namespace seqan3::execution
{

template <typename sender_t, typename receiver_t>
using connect_result_t = std::invoke_result_t<decltype(connect), sender_t, receiver_t>;
}  // namespace seqan3::execution
