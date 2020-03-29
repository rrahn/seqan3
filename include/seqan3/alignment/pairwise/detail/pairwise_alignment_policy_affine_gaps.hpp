// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_policy_affine_gaps.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/pairwise/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>

namespace seqan3::detail
{

template <typename alignment_config_t, bool with_zero_cutoff>
class pairwise_alignment_policy_affine_gaps
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using score_type = typename traits_type::score_t;
    using affine_cell_proxy_t = affine_cell_proxy<std::tuple<score_type, score_type, score_type>>;

protected:
    //!\brief The score for a gap extension.
    score_type gap_extension_score{};
    //!\brief The score for a gap opening including the gap extension.
    score_type gap_open_score{};
    //!\brief Flag indicating whether the first row is zero initialised.
    bool zero_initialise_row{};
    //!\brief Flag indicating whether the first column is zero initialised.
    bool zero_initialise_column{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_policy_affine_gaps() = default; //!< Defaulted.
    pairwise_alignment_policy_affine_gaps(pairwise_alignment_policy_affine_gaps const &) = default; //!< Defaulted.
    pairwise_alignment_policy_affine_gaps(pairwise_alignment_policy_affine_gaps &&) = default; //!< Defaulted.
    pairwise_alignment_policy_affine_gaps & operator=(pairwise_alignment_policy_affine_gaps const &) = default; //!< Defaulted.
    pairwise_alignment_policy_affine_gaps & operator=(pairwise_alignment_policy_affine_gaps &&) = default; //!< Defaulted.
    ~pairwise_alignment_policy_affine_gaps() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param config The alignment configuration.
     *
     * \details
     *
     * Initialises the gap open score and gap extension score for this policy.
     */
    explicit pairwise_alignment_policy_affine_gaps(alignment_config_t const & config)
    {
        // Get the gap scheme from the config and choose -1 and -10 as default.
        gap_scheme my_gap{gap_score{-1}, seqan3::gap_open_score{-10}};
        auto && tmp_scheme = config.template value_or<align_cfg::gap>(my_gap);
        gap_extension_score = static_cast<score_type>(tmp_scheme.get_gap_score());
        gap_open_score = static_cast<score_type>(tmp_scheme.get_gap_open_score()) + gap_extension_score;

        auto align_ends_cfg = config.template value_or<align_cfg::aligned_ends>(free_ends_none);
        zero_initialise_row = align_ends_cfg[0];
        zero_initialise_column = align_ends_cfg[2];
    }
    //!\}

    /*!\brief Initialises the first cell of the alignment matrix in the top left corner of the matrix.
     * \tparam affine_cell_t The type of the affine cell.
     *
     * \param[in] diagonal_score The previous diagonal score, which corresponds to \f$M[i - 1, j - 1]\f$.
     * \param[in] previous_cell The predecessor cell corresponding to the values \f$V[i - 1, j]\f$ and \f$H[i, j -1]\f$.
     *
     * \returns The computed cell.
     *
     * \details
     *
     * Computes the current cell according to following recursion formula:
     * \f[
     * H[i, j] = \max
     * \begin{cases}
     *  M[i, j - 1] + g_o\\
     *  H[i, j - 1] + g_e
     * \end{cases}
     *
     * V[i, j] = \max
     * \begin{cases}
     *  M[i - 1, j] + g_o\\
     *  V[i - 1, j] + g_e
     * \end{cases}
     *
     * M[i, j] = \max
     * \begin{cases}
     *  M[i - 1, j - 1] + \delta(s1_i, s2_j)\\
     *  V[i, j] \\
     *  H[i, j]
     * \end{cases}
     * \f]
     */
    template <typename affine_cell_t>
    auto compute_inner_cell(score_type diagonal_score,
                            affine_cell_t previous_cell,
                            score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();

        // std::cout << diagonal_score << ": " <<   tmp << " ";
        diagonal_score = (diagonal_score < vertical_score) ? vertical_score : diagonal_score;
        diagonal_score = (diagonal_score < horizontal_score) ? horizontal_score : diagonal_score;

        // Check if zero cutoff is enabled so that score cannot fall below 0.
        if constexpr (with_zero_cutoff)
            diagonal_score = (diagonal_score < 0) ? 0 : diagonal_score;

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        // store the vertical_score and horizontal_score value in the next path
        vertical_score = (vertical_score < tmp) ? tmp : vertical_score;
        horizontal_score = (horizontal_score < tmp) ? tmp : horizontal_score;

        return affine_cell_proxy_t{diagonal_score, horizontal_score, vertical_score};
    }

    /*!\brief Initialises the first cell of the alignment matrix in the top left corner of the matrix.
     *
     * \returns The computed cell.
     *
     * \details
     *
     * Initialises the cell at the origin of the alignment matrix (top left corner of the matrix). The optimal score is
     * initialised to 0, while the value of the horizontal and vertical matrix are initialised as:
     * \f$V[0, 0] = H[0, 0] = g_o\f$.
     */
    affine_cell_proxy_t initialise_origin_cell() const noexcept
    {
        return {0, (zero_initialise_row) ? 0 : gap_open_score, (zero_initialise_column) ? 0 : gap_open_score};
    }

    /*!\brief Initialises a cell of the first alignment matrix column.
     * \tparam affine_cell_t The type of the affine cell.
     *
     * \param previous_cell The predecessor cell on the same column.
     *
     * \returns The computed cell.
     *
     * \details
     *
     * Initialises a cell of the first alignment matrix column. The optimal score is the same as the vertical score
     * which is equal to \f$V[i, 0] = M[i, 0] = g_o + g_e * i\f$. The horizontal score is initialised to
     * \f$H[i, 0] = V[i, 0] + g_o\f$ to prohibit extending a gap in the horizontal matrix from \f$H[i, 0]\f$.
     */
    template <typename affine_cell_t>
    affine_cell_proxy_t initialise_first_column_cell(affine_cell_t cell) const noexcept
    {
        score_type vertical_score = (zero_initialise_column) ? 0 : cell.vertical_score() + gap_extension_score;
        return {cell.vertical_score(), cell.vertical_score() + gap_open_score, vertical_score};
    }

    /*!\brief Initialises the first cell of a alignment matrix column.
     * \tparam affine_cell_t The type of the affine cell.
     *
     * \param cell The predecessor cell on the same row.
     *
     * \returns The computed cell.
     *
     * \details
     *
     * Initialises the first cell of a alignment matrix column. The optimal score is the same as the horizontal score
     * which is equal to \f$H[0, j] = M[0, j] = g_o + g_e * j\f$. The vertical score is initialised to
     * \f$H[0, j] + g_o\f$ to prohibit extending a gap in the vertical matrix from \f$V[0, j]\f$.
     */
    template <typename affine_cell_t>
    affine_cell_proxy_t initialise_first_row_cell(affine_cell_t cell) const noexcept
    {
        score_type horizontal_score = (zero_initialise_row) ? 0 : cell.horizontal_score() + gap_extension_score;
        return {cell.horizontal_score(),
                horizontal_score,
                cell.horizontal_score() + gap_open_score};
    }
};
} // namespace seqan3::detail
