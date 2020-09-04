// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_recursion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the alignment algorithm using affine gap costs.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration.
 *
 * \details
 *
 * Implements the functions to initialise and compute the alignment matrix using the recursion formula for affine gaps.
 * Other policies can inherit from this policy and overload the recursion functions, e.g. to change the
 * initialisation of the alignment matrix.
 *
 * \note For more information, please refer to the original article for the alignment with affine gap cost function:
 *       GOTOH, Osamu. An improved algorithm for matching biological sequences.
 *       Journal of molecular biology, 1982, 162. Jg., Nr. 3, S. 705-708.
 */
template <typename alignment_configuration_t>
class policy_affine_gap_recursion
{
protected:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;
    //!\brief The configured trace type.
    using trace_type = typename traits_type::trace_type;
    //!\brief The internal tuple storing the scores of an affine cell.
    using affine_score_tuple_t = std::tuple<score_type, score_type, score_type>;
    //!\brief The internal tuple storing the trace directions of an affine cell.
    using affine_trace_tuple_t = std::tuple<trace_type, trace_type, trace_type>;
    //!\brief The affine cell type returned by the functions.
    using affine_cell_type = std::conditional_t<traits_type::requires_trace_information,
                                                affine_cell_proxy<std::pair<affine_score_tuple_t,
                                                                            affine_trace_tuple_t>>,
                                                affine_cell_proxy<affine_score_tuple_t>>;

    //!\brief The score for a gap extension.
    score_type gap_extension_score{};
    //!\brief The score for a gap opening including the gap extension.
    score_type gap_open_score{};

    //!\brief Initialisation state of the first row of the alignment.
    bool first_row_is_free{};
    //!\brief Initialisation state of the first column of the alignment.
    bool first_column_is_free{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_recursion() = default; //!< Defaulted.
    policy_affine_gap_recursion(policy_affine_gap_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_recursion(policy_affine_gap_recursion &&) = default; //!< Defaulted.
    policy_affine_gap_recursion & operator=(policy_affine_gap_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_recursion & operator=(policy_affine_gap_recursion &&) = default; //!< Defaulted.
    ~policy_affine_gap_recursion() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration.
     *
     * \details
     *
     * Initialises the gap open score and gap extension score for this policy.
     * If no gap cost model was provided by the user the default gap costs `-10` and `-1` are set for the gap open score
     * and the gap extension score respectively.
     */
    explicit policy_affine_gap_recursion(alignment_configuration_t const & config)
    {
        // Get the gap scheme from the config or choose -1 and -10 as default.
        auto const & selected_gap_scheme = config.get_or(align_cfg::gap{gap_scheme{seqan3::gap_score{-1},
                                                                                   seqan3::gap_open_score{-10}}}).value;

        if constexpr (simd::simd_concept<score_type>)
        {
            gap_extension_score = simd::fill<score_type>(selected_gap_scheme.get_gap_score());
            gap_open_score = simd::fill<score_type>(selected_gap_scheme.get_gap_open_score()) + gap_extension_score;
        }
        else
        {
            gap_extension_score = static_cast<score_type>(selected_gap_scheme.get_gap_score());
            gap_open_score = static_cast<score_type>(selected_gap_scheme.get_gap_open_score()) + gap_extension_score;
        }

        auto method_global_config = config.get_or(align_cfg::method_global{});
        first_row_is_free = method_global_config.free_end_gaps_sequence1_leading;
        first_column_is_free = method_global_config.free_end_gaps_sequence2_leading;
    }
    //!\}

    /*!\brief Computes an inner cell of the alignment matrix.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] diagonal_score The previous diagonal score, which corresponds to \f$M[i - 1, j - 1]\f$.
     * \param[in] previous_cell The predecessor cell corresponding to the values \f$V[i - 1, j]\f$ and \f$H[i, j -1]\f$.
     * \param[in] sequence_score The score obtained from the scoring scheme for the current cell (\f$ \delta\f$).
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Computes the current cell according to following recursion formula:
     * * \f$ H[i, j] = \max \{M[i, j - 1] + g_o, H[i, j - 1] + g_e\}\f$
     * * \f$ V[i, j] = \max \{M[i - 1, j] + g_o, V[i - 1, j] + g_e\}\f$
     * * \f$ M[i, j] = \max \{M[i - 1, j - 1] + \delta, H[i, j], V[i, j]\}\f$
     */
    template <affine_cell_proxy_instance affine_cell_t>
    affine_cell_type compute_inner_cell(score_type diagonal_score,
                                        affine_cell_t previous_cell,
                                        score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();

        diagonal_score = (diagonal_score < vertical_score) ? vertical_score : diagonal_score;
        diagonal_score = (diagonal_score < horizontal_score) ? horizontal_score : diagonal_score;

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        // store the vertical_score and horizontal_score value in the next path
        vertical_score = (vertical_score < tmp) ? tmp : vertical_score;
        horizontal_score = (horizontal_score < tmp) ? tmp : horizontal_score;

        return {diagonal_score, horizontal_score, vertical_score};
    }

    //!\overload
    template <affine_cell_proxy_with_trace_instance affine_cell_t>
    affine_cell_type compute_inner_cell(score_type diagonal_score,
                                        affine_cell_t previous_cell,
                                        score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();
        trace_directions best_trace = trace_directions::diagonal;

        diagonal_score = (diagonal_score < vertical_score)
                       ? (best_trace = previous_cell.vertical_trace(), vertical_score)
                       : (best_trace |= previous_cell.vertical_trace() ,diagonal_score);
        diagonal_score = (diagonal_score < horizontal_score)
                       ? (best_trace = previous_cell.horizontal_trace(), horizontal_score)
                       : (best_trace |= previous_cell.horizontal_trace(), diagonal_score);

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        // store the vertical_score and horizontal_score value in the next path
        trace_directions next_vertical_trace = trace_directions::up;
        trace_directions next_horizontal_trace = trace_directions::left;

        vertical_score = (vertical_score < tmp)
                       ? (next_vertical_trace = trace_directions::up_open, tmp)
                       : vertical_score;
        horizontal_score = (horizontal_score < tmp)
                         ? (next_horizontal_trace = trace_directions::left_open, tmp)
                         : horizontal_score;

        return {{diagonal_score, horizontal_score, vertical_score},
                {best_trace, next_horizontal_trace, next_vertical_trace}};
    }

    /*!\brief Initialises the first cell of the alignment matrix in the top left corner of the matrix.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises the cell at the origin of the alignment matrix (top left corner of the matrix). The optimal score is
     * initialised to 0, while the value of the horizontal and vertical matrix are initialised as:
     * \f$V[0, 0] = H[0, 0] = g_o\f$.
     */
    affine_cell_type initialise_origin_cell() const noexcept
        requires affine_cell_proxy_instance<affine_cell_type>
    {
        return {score_type{},
                first_row_is_free ? score_type{} : gap_open_score,
                first_column_is_free ? score_type{} : gap_open_score};
    }

    //!\overload
    affine_cell_type initialise_origin_cell() const noexcept
        requires affine_cell_proxy_with_trace_instance<affine_cell_type>
    {
        return {{score_type{},
                 first_row_is_free ? score_type{} : gap_open_score,
                 first_column_is_free ? score_type{} : gap_open_score},
                {trace_directions::none,
                 first_row_is_free ? trace_directions::none : trace_directions::left_open,
                 first_column_is_free ? trace_directions::none : trace_directions::up_open}};
    }

    /*!\brief Initialises a cell of the first alignment matrix column.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] previous_cell The predecessor cell on the same column \f$M[i-1, 0]\f$.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises a cell of the first alignment matrix column. The optimal score is the same as the vertical score
     * which is equal to \f$V[i, 0] = M[i, 0] = g_o + g_e * i\f$. The horizontal score is initialised to
     * \f$H[i, 0] = V[i, 0] + g_o\f$ to prohibit extending a gap in the horizontal matrix from \f$H[i, 0]\f$.
     */
    template <affine_cell_proxy_instance affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_vertical = previous_cell.vertical_score() + gap_extension_score;
        return {previous_cell.vertical_score(),
                previous_cell.vertical_score() + gap_open_score,
                first_column_is_free ? previous_cell.vertical_score() : new_vertical};
    }

    //!\overload
    template <affine_cell_proxy_with_trace_instance affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_vertical = previous_cell.vertical_score() + gap_extension_score;
        return {{previous_cell.vertical_score(),
                 previous_cell.vertical_score() + gap_open_score,
                 first_column_is_free ? previous_cell.vertical_score() : new_vertical},
                {previous_cell.vertical_trace(),
                 trace_directions::left_open,
                 first_column_is_free ? trace_directions::none : trace_directions::up}};
    }

    /*!\brief Initialises the first cell of a alignment matrix column.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] previous_cell The predecessor cell on the same row \f$M[0, j-1]\f$.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises the first cell of a alignment matrix column. The optimal score is the same as the horizontal score
     * which is equal to \f$H[0, j] = M[0, j] = g_o + g_e * j\f$. The vertical score is initialised to
     * \f$V[0,j] = H[0, j] + g_o\f$ to prohibit extending a gap in the vertical matrix from \f$V[0, j]\f$.
     */
    template <affine_cell_proxy_instance affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_horizontal_score = previous_cell.horizontal_score() + gap_extension_score;
        return {previous_cell.horizontal_score(),
                first_row_is_free ? previous_cell.horizontal_score() : new_horizontal_score,
                previous_cell.horizontal_score() + gap_open_score};
    }

    //!\overload
    template <affine_cell_proxy_with_trace_instance affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_horizontal_score = previous_cell.horizontal_score() + gap_extension_score;
        return {{previous_cell.horizontal_score(),
                 first_row_is_free ? previous_cell.horizontal_score() : new_horizontal_score,
                 previous_cell.horizontal_score() + gap_open_score},
                {previous_cell.horizontal_trace(),
                 first_row_is_free ? trace_directions::none : trace_directions::left,
                 trace_directions::up_open}};
    }

    /*!\brief Returns the lowest viable score.
     *
     * \details
     *
     * In some versions of the algorithms a value representing minus infinity is needed. Since the data type is an
     * signed integral there is no infinity but only the lowest possible value that can be represented by the score
     * type. In order to avoid unnecessary if conditions to protect against signed integer underflow the lowest viable
     * score is computed. Subtracting a gap penalty from this will still result in a valid score which represents
     * minus infinity.
     */
    score_type lowest_viable_score() const noexcept
    {
        assert(gap_open_score <= 0 && gap_extension_score <= 0);
        return std::numeric_limits<score_type>::lowest() - (gap_open_score + gap_extension_score);
    }
};
} // namespace seqan3::detail
