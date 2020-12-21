// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::operation_state concept.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/execution/cpo/start.hpp>

namespace seqan3::execution
{

template<class operation_t>
SEQAN3_CONCEPT operation_state = destructible<operation_t> &&
                                 is_object_v<operation_t> &&
                                 requires (operation_t & o)
                                 {
                                     { execution::start(o) } noexcept;
                                 };
}  // namespace seqan3::execution
