// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::sender_traits.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

template <template <template <typename ...> typename, template <typename ...> typename> typename>
struct has_value_types;

template <template <template <typename ...> typename> typename>
struct has_error_types;

template<typename sender_t>
SEQAN3_CONCEPT has_sender_types = requires
{
    typename has_value_types <sender_t::template value_types>;
    typename has_error_types <sender_t::template error_types>;
    typename bool_constant<sender_t::sends_done>;
};

} // namespace seqan3::detail

namespace seqan3::execution
{

template <detail::has_sender_types sender_t>
struct sender_traits
{
    template <template <typename ...> typename Tuple, template <typename ...> typename Variant>
    using value_types = typename sender_t::template value_types<Tuple, Variant>;

    template <template <typename ...> typename Variant>
    using error_types = typename sender_t::template error_types<Variant>;

    static constexpr bool sends_done = sender_t::sends_done;
};

template <typename sender_t>
struct sender_traits<sender_t>
{
    template <template <typename ...> typename Tuple, template <typename ...> typename Variant>
    using value_types = Variant<Tuple<>>;

    template <template <typename ...> typename Variant>
    using error_types = Variant<std::exception_ptr>;

    static constexpr bool sends_done = true;
};
}  // namespace seqan3::execution
