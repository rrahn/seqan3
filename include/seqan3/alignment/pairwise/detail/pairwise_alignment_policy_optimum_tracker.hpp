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

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

template <bool activate_last_row_tracking, bool activate_last_column_tracking, bool activate_all_tracking>
struct optimum_tracker_traits
{
    static constexpr bool track_last_row_cell = activate_last_row_tracking;
    static constexpr bool track_last_column_cell = activate_last_column_tracking;
    static constexpr bool track_every_cell = activate_all_tracking;
};

template <typename alignment_config_t, typename tracker_traits_t = optimum_tracker_traits<false, false, false>>
class pairwise_alignment_policy_optimum_tracker
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using score_type = typename traits_type::score_t;

    score_type optimal_score{std::numeric_limits<score_type>::lowest()};
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

    pairwise_alignment_policy_optimum_tracker(alignment_config_t const & /*config*/)
    {}
    //!\}

    template <typename cell_t>
    decltype(auto) track_cell(cell_t && cell) noexcept
    {
        track_cell<tracker_traits_t::track_every_cell>(cell);
        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_last_row_cell(cell_t && cell) noexcept
    {
        track_cell<tracker_traits_t::track_last_row_cell>(cell);
        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_last_column_cell(cell_t && cell) noexcept
    {
        track_cell<tracker_traits_t::track_last_column_cell>(cell);
        return std::forward<cell_t>(cell);
    }

    template <typename cell_t>
    decltype(auto) track_final_cell(cell_t && cell) noexcept
    {
        track_cell<true>(cell);
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

    template <bool check_optimum, typename cell_t>
    void track_cell(cell_t && cell) noexcept
    {
        if constexpr (check_optimum)
            optimal_score = (cell.optimal_score() >= optimal_score) ? cell.optimal_score() : optimal_score;
    }
};
} // namespace seqan3::detail
