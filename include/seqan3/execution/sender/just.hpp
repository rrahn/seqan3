// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::execution::just.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <seqan3/std/type_traits>

#include <seqan3/execution/cpo/set_value.hpp>
#include <seqan3/execution/cpo/set_error.hpp>

namespace seqan3::execution::detail
{

template <typename ...send_values_t>
class sender_just
{
    using values_t = std::tuple<send_values_t...>;
    // We cannot simply forward references.
    values_t stored_values;

    template <typename receiver_t>
    class just_operation_state
    {
    private:
        values_t stored_values;
        std::remove_cvref_t<receiver_t> output_receiver;

    public:

        just_operation_state(values_t values, receiver_t receiver) :
            stored_values{std::move(values)},
            output_receiver{std::move(receiver)}
        {}

        void start() noexcept
        {
            try
            {
                std::apply(
                    [this] (auto && ...values)
                    {
                        execution::set_value(std::move(output_receiver), std::move(values)...);
                    }, std::move(stored_values));
            }
            catch (...)
            {
                execution::set_error(std::move(output_receiver), std::current_exception());
            }
        }
    };

public:

    template <template <typename...> class Tuple, template <typename...> class Variant>
    using value_types = Variant<Tuple<send_values_t...>>;

    template <template <typename...> class Variant>
    using error_types = Variant<std::exception_ptr>;

    static constexpr bool sends_done = false;

    sender_just() = default;

    template <typename ...output_values_t>
    sender_just(output_values_t && ...values) : stored_values{std::move(values)...}
    {}

    template <typename receiver_t>
    auto connect(receiver_t && output_receiver) && noexcept(std::is_nothrow_move_constructible_v<values_t>)
    {
        return just_operation_state{std::move(stored_values), std::forward<receiver_t>(output_receiver)};
    }

    template <typename receiver_t>
    auto connect(receiver_t && output_receiver) const & noexcept(std::is_nothrow_copy_constructible_v<values_t>)
    {
        return just_operation_state{stored_values, std::forward<receiver_t>(output_receiver)};
    }
};

} // namespace seqan3::execution::detail

namespace seqan3::execution
{
inline constexpr auto just = [] (auto && ...values)
{
    return detail::sender_just<std::remove_cvref_t<decltype(values)>...>{std::forward<decltype(values)>(values)...};
};
}  // namespace seqan3::execution
