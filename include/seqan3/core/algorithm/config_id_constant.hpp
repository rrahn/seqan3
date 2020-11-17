// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::config_id_constant.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <type_traits>

namespace seqan3::detail
{

template <typename config_id_t>
//!\cond
    requires std::is_enum_v<config_id_t>
//!\endcond
class config_id_constant
{
private:
    uint64_t value{};

public:
    //!\brief The type of the configuration id.
    using type = config_id_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr config_id_constant() = default; //!< Defaulted.
    constexpr config_id_constant(config_id_constant const &) = default; //!< Defaulted.
    constexpr config_id_constant(config_id_constant &&) = default; //!< Defaulted.
    constexpr config_id_constant & operator=(config_id_constant const &) = default; //!< Defaulted.
    constexpr config_id_constant & operator=(config_id_constant &&) = default; //!< Defaulted.
    ~config_id_constant() = default; //!< Defaulted.

    template <typename ...config_id_type>
        requires (std::same_as<config_id_type, config_id_t> && ...)
    constexpr config_id_constant(config_id_type const ...ids) noexcept
    {
        value = ((1ull << to_integral(ids)) | ...);
    }
    //!\}

    /*!\brief Tests whether the id is contained in the config_id_constant.
     *
     * \param[in] id The id to test for.
     *
     * \returns `true` if the id is contained, otherwise `false`.
     *
     * \details
     *
     * Checks if the given id i
     *
     */
    constexpr bool contains(config_id_t const id) const noexcept
    {
        assert((static_cast<uint64_t>(id) + 1) < (sizeof(uint64_t) * 8));

        return value & (1ull << to_integral(id));
    }
private:
    constexpr uint64_t to_integral(config_id_t const id) const noexcept
    {
        return static_cast<uint64_t>(id) + 1;
    }
};

}  // namespace seqan3::detail
