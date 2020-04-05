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
#include <seqan3/alignment/pairwise/detail/alignment_coordinate_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/type_traits/basic.hpp>

namespace seqan3::detail
{

template <typename alignment_config_t, bool track_every_cell = false>
class pairwise_alignment_policy_optimum_tracker
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using score_type = typename traits_type::score_t;

    using coordinate_type = typename alignment_coordinate_matrix<true>::alignment_coordinate;

    static constexpr score_type infinity = [] () constexpr
    {
        if constexpr (simd_concept<score_type>)
            return simd::fill<score_type>(std::numeric_limits<typename simd_traits<score_type>::scalar_type>::lowest());
        else
            return std::numeric_limits<score_type>::lowest();
    }();

    score_type optimal_score{infinity};
    coordinate_type m_coordinate{std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::lowest()};
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

    template <typename cell_t, typename coordinate_t>
    decltype(auto) track_cell(cell_t && cell, [[maybe_unused]] coordinate_t && coordinate) noexcept
    {
        if constexpr (track_every_cell)
            exchange_max(cell, std::forward<coordinate_t>(coordinate));

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t, typename coordinate_t>
    decltype(auto) track_last_row_cell(cell_t && cell, [[maybe_unused]] coordinate_t && coordinate) noexcept
    {
        if (track_last_row)
            exchange_max(cell, std::forward<coordinate_t>(coordinate));

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t, typename coordinate_t>
    decltype(auto) track_last_column_cell(cell_t && cell, [[maybe_unused]] coordinate_t && coordinate) noexcept
    {
        if (track_last_column)
            exchange_max(cell, std::forward<coordinate_t>(coordinate));

        return std::forward<cell_t>(cell);
    }

    template <typename cell_t, typename coordinate_t>
    decltype(auto) track_final_cell(cell_t && cell, [[maybe_unused]] coordinate_t && coordinate) noexcept
    {
        exchange_max(cell, std::forward<coordinate_t>(coordinate));
        return std::forward<cell_t>(cell);
    }

    auto const & optimum_score() const noexcept
    {
        return optimal_score;
    }

    coordinate_type const & optimum_coordinate() const noexcept
    {
        return m_coordinate;
    }

    void reset_tracker() noexcept
    {
        optimal_score = infinity;
        m_coordinate = coordinate_type{std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::lowest()};
    }

private:

    template <alignment_score_cell cell_t, typename coordinate_t>
    void exchange_max(cell_t && cell, coordinate_t && coordinate) noexcept
    {
        if constexpr (decays_to_ignore_v<std::remove_reference_t<coordinate_t>>)
            optimal_score = (cell.optimal_score() >= optimal_score) ? cell.optimal_score() : optimal_score;
        else // update the coordinate if a new max was found and coordinate is enabled.
            optimal_score = (cell.optimal_score() >= optimal_score)
                          ? (m_coordinate = coordinate, cell.optimal_score())
                          : optimal_score;
    }

    template <alignment_score_trace_cell cell_t, typename coordinate_t>
    void exchange_max(cell_t && cell, coordinate_t && coordinate) noexcept
    {
        using std::get;

        if constexpr (decays_to_ignore_v<std::remove_reference_t<coordinate_t>>)
            optimal_score = (get<0>(cell).optimal_score() >= optimal_score)
                          ? get<0>(cell).optimal_score()
                          : optimal_score;
        else // update the coordinate if a new max was found and coordinate is enabled.
            optimal_score = (get<0>(cell).optimal_score() >= optimal_score)
                          ? (m_coordinate = coordinate, get<0>(cell).optimal_score())
                          : optimal_score;
    }
};
} // namespace seqan3::detail
