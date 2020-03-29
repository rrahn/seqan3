// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_policy_optimum_tracker.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>

namespace seqan3::detail
{

template <typename alignment_config_t, bool track_every_cell = false>
class pairwise_alignment_policy_optimum_tracker
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using score_type = typename traits_type::score_t;

    score_type optimal_score{std::numeric_limits<score_type>::lowest()};
    bool track_last_row{};
    bool track_last_column{};
protected:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    pairwise_alignment_policy_optimum_tracker() = default;
    //!\brief Defaulted.
    pairwise_alignment_policy_optimum_tracker(pairwise_alignment_policy_optimum_tracker const &) = default;
    //!\brief Defaulted.
    pairwise_alignment_policy_optimum_tracker(pairwise_alignment_policy_optimum_tracker &&) = default;
    //!\brief Defaulted.
    pairwise_alignment_policy_optimum_tracker & operator=(pairwise_alignment_policy_optimum_tracker const &) = default;
    //!\brief Defaulted.
    pairwise_alignment_policy_optimum_tracker & operator=(pairwise_alignment_policy_optimum_tracker &&) = default;
    //!\brief Defaulted.
    ~pairwise_alignment_policy_optimum_tracker() = default;

    pairwise_alignment_policy_optimum_tracker(alignment_config_t const & config)
    {
        auto align_ends_cfg = config.template value_or<align_cfg::aligned_ends>(free_ends_none);
        track_last_row = align_ends_cfg[1];
        track_last_column = align_ends_cfg[3];
    }
    //!\}

    template <typename cell_t>
    decltype(auto) track_cell(cell_t && cell) noexcept
    {
        if constexpr (track_every_cell)
            exchange_max(cell);

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_last_row_cell(cell_t && cell) noexcept
    {
        if (track_last_row)
            exchange_max(cell);

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_last_column_cell(cell_t && cell) noexcept
    {
        if (track_last_column)
            exchange_max(cell);

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_final_cell(cell_t && cell) noexcept
    {
        exchange_max(cell);
        return std::forward<cell_t>(cell);
    }

    score_type optimum() const noexcept
    {
        return optimal_score;
    }

    void reset_tracker() noexcept
    {
        optimal_score = std::numeric_limits<score_type>::lowest();
    }

private:

    template <typename cell_t>
    void exchange_max(cell_t && cell) noexcept
    {
        optimal_score = (cell.optimal_score() >= optimal_score) ? cell.optimal_score() : optimal_score;
    }
};
} // namespace seqan3::detail
