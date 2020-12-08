// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_optimum_tracker.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <limits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief A function object that compares and possibly updates the alignment optimum with the current cell.
 * \ingroup pairwise_alignment
 *
 * \details
 *
 * Updates the current alignment optimum with the new score and the respective coordinate if the new score
 * compares greater or equal to the score of the current optimum.
 *
 * \see seqan3::detail::max_score_banded_updater
 */
struct max_score_updater
{
    /*!\brief Compares and updates the optimal score-coordinate pair.
     * \tparam score_t The type of the score to track; must model std::totally_ordered and
     *                 std::assignable_from `const & score_t`.
     * \tparam coordinate_t The type of the coordinate to track; must model std::assignable_from `const & coordinate_t`.
     *
     * \param[in,out] optimal_score The optimal score to update.
     * \param[in,out] optimal_coordinate The optimal coordinate to update.
     * \param[in] current_score The score of the current cell.
     * \param[in] current_coordinate The coordinate of the current cell.
     *
     * \details
     *
     * Compares the current_score with the optimal score and updates the optimal score and coordinate if the current
     * one is the new optimum. Otherwise, keeps the old optimum.
     */
    template <typename score_t, typename coordinate_t>
    //!\cond
        requires (std::totally_ordered<score_t> && std::assignable_from<score_t &, score_t const &> &&
                  std::assignable_from<coordinate_t &, coordinate_t const &>)
    //!\endcond
    void operator()(score_t & optimal_score,
                    coordinate_t & optimal_coordinate,
                    score_t current_score,
                    coordinate_t current_coordinate) const noexcept
    {
        bool const is_better_score = current_score > optimal_score;
        optimal_score = (is_better_score) ? std::move(current_score) : optimal_score;
        optimal_coordinate = (is_better_score) ? std::move(current_coordinate) : optimal_coordinate;
    }
};

/*!\brief A function object that compares and possibly updates the alignment optimum with the current cell.
 * \ingroup pairwise_alignment
 *
 * \details
 *
 * This special function object is used for global alignments including free end-gaps.
 * It updates the current alignment optimum with the new score and the respective coordinate if the new score
 * is greater or equal to the score of the current optimum and the cell is either in the target row or column.
 * For the banded alignment this additional check is needed because the band might not reach the last row of the matrix
 * for all columns.
 * Thus, the last cell of a column will only be checked if the corresponding coordinate points to the last row of the
 * matrix.
 *
 * ### Example
 *
 * Consider the following banded matrix:
 *
 * ```
 *    0         1
 *    012345678901
 *    ____________
 * 0 |       \    |
 * 1 |        \   |
 * 2 |\        \  |
 * 3 | \        \ |
 * 4 |  \        *|
 * 5 |   \       *|
 * 6 |    ********|
 *    ––––––––––––
 * ```
 *
 * When computing the columns from index 0 to 3 the last row of the matrix is not covered by the band.
 * But the algorithm that computes the column does not know if the last computed cell of it is in fact also the last
 * cell of the matrix (row index 6). Thus, by comparing the cell coordinate, which corresponds to the absolute
 * matrix coordinate, with the initialised target row and column index of this function object, it can be guaranteed
 * that only the scores of the cells are tracked that are in the last row or column (cells marked with an `*`
 * in the picture above).
 *
 * \see seqan3::detail::max_score_updater
 */
struct max_score_banded_updater
{
private:
    //!\brief The row index indicating the last row of the alignment matrix.
    size_t target_row_index{};
    //!\brief The column index indicating the last column of the alignment matrix.
    size_t target_col_index{};

public:
    /*!\brief Compares and updates the optimal score-coordinate pair.
     * \tparam score_t The type of the score to track; must model std::totally_ordered and
     *                 std::assignable_from `const & score_t`.
     * \tparam coordinate_t The type of the coordinate to track; must model std::assignable_from `const & coordinate_t`.
     *
     * \param[in,out] optimal_score The optimal score to update.
     * \param[in,out] optimal_coordinate The optimal coordinate to update.
     * \param[in] current_score The score of the current cell.
     * \param[in] current_coordinate The coordinate of the current cell.
     *
     * \details
     *
     * First, checks if the coordinate of the current alignment cell is either in the target row or column, i.e.
     * the last row/column of the alignment matrix. If this is not the case the score won't be updated.
     * Otherwise, compares the current_score with the optimal score and updates the optimal score and coordinate
     * accordingly.
     */
    template <typename score_t, typename coordinate_t>
    //!\cond
        requires (std::totally_ordered<score_t> && std::assignable_from<score_t &, score_t const &> &&
                  std::assignable_from<coordinate_t &, coordinate_t const &>)
    //!\endcond
    void operator()(score_t & optimal_score,
                    coordinate_t & optimal_coordinate,
                    score_t current_score,
                    coordinate_t current_coordinate) const noexcept
    {
        bool const is_better_score = (target_row_index == current_coordinate.row ||
                                      target_col_index == current_coordinate.col) &&
                                      (current_score >= optimal_score);
        optimal_score = (is_better_score) ? std::move(current_score) : optimal_score;
        optimal_coordinate = (is_better_score) ? std::move(current_coordinate) : optimal_coordinate;
    }

    /*!\brief Sets the target index for the last row and column of the matrix.
     * \param row_index The target index of the last row.
     * \param col_index The target index of the last column.
     */
    void set_target_indices(row_index_type<size_t> row_index, column_index_type<size_t> col_index) noexcept
    {
        target_row_index = row_index.get();
        target_col_index = col_index.get();
    }
};

/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup pairwise_alignment
 *
 * \tparam score_t The type of the optimal score to track; must model seqan3::arithmetic or
 *                 seqan3::simd::simd_concept.
 * \tparam matrix_coordinate_t The type of the matrix coordinate to track; must model
 *                             seqan3::detail::template_specialisation_of seqan3::detail::matrix_index.
 * \tparam optimum_updater_t The type of the optimum update operation, which compares and updates the alignment optimum
 *                           with the current cell; must model std::semiregular.
 * \details
 *
 * Implements the interface to track the alignment optimum. It updates the currently stored optimum using the
 * optimum update operation. The optimum update operation is stored inside of the class and can have a state.
 * The optimum updater must be invokable with a reference to the optimal score and coordinate and the
 * score and coordinate of the current cell.
 * Special methods are offered to track any cell (for example when computing the local alignment), the last cell of a
 * column or a row (for example when using free-end gaps), or the final cell of the entire matrix (for example in the
 * standard global alignment).
 * The optimum needs to be reset in between alignment computations in order to ensure that the correct result is
 * tracked.
 */
template <typename score_t, typename matrix_coordinate_t, std::semiregular optimum_updater_t>
//!\cond
    requires arithmetic<score_t> || simd_concept<score_t> &&
             is_type_specialisation_of_v<matrix_coordinate_t, matrix_index> &&
             std::invocable<optimum_updater_t, score_t &, matrix_coordinate_t &, score_t, matrix_coordinate_t>
//!\endcond
class policy_optimum_tracker
{
public:
    // ----------------------------------------------------------------------------
    // type definitions
    // ----------------------------------------------------------------------------

    //!\brief The configured score type.
    using score_type = score_t;
    //!\brief The matrix coordinate type that is used to locate a cell inside of the alignment matrix.
    using matrix_coordinate_type = matrix_coordinate_t;

    // ----------------------------------------------------------------------------
    // non-static members
    // ----------------------------------------------------------------------------

    //!\brief The tracked score of the global optimum.
    score_type optimal_score{};
    //!\brief The matrix coordinate of the tracked optimum.
    matrix_coordinate_type optimal_coordinate{};
    //!\brief The function object to compare and exchange the optimum.
    optimum_updater_t compare_and_set_optimum{};

    //!\brief Whether every cell of the alignment matrix shall be tracked.
    bool test_every_cell{false};
    //!\brief Whether cells of the last row shall be tracked.
    bool test_last_row_cell{false};
    //!\brief Whether cells of the last column shall be tracked.
    bool test_last_column_cell{false};

    // ----------------------------------------------------------------------------
    // member functions
    // ----------------------------------------------------------------------------

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_optimum_tracker() = default; //!< Defaulted.
    policy_optimum_tracker(policy_optimum_tracker const &) = default; //!< Defaulted.
    policy_optimum_tracker(policy_optimum_tracker &&) = default; //!< Defaulted.
    policy_optimum_tracker & operator=(policy_optimum_tracker const &) = default; //!< Defaulted.
    policy_optimum_tracker & operator=(policy_optimum_tracker &&) = default; //!< Defaulted.
    ~policy_optimum_tracker() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration (not used in this context).
     *
     * \details
     *
     * Reads the state of seqan3::align_cfg::method_global and enables the tracking of the last row or column if
     * requested. Otherwise, only the last cell will be tracked.
     */
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    policy_optimum_tracker(alignment_configuration_t const & config)
    {
        auto method_global_config = config.get_or(align_cfg::method_global{});
        constexpr bool track_all_cells = config.template exists<align_cfg::method_local>();
        test_last_row_cell = method_global_config.free_end_gaps_sequence1_trailing || track_all_cells;
        test_last_column_cell = method_global_config.free_end_gaps_sequence2_trailing || track_all_cells;
        test_every_cell = track_all_cells;
    }
    //!\}

    /*!\brief Tracks any cell within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `best_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a local alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_every_cell)
            invoke_comparator(cell, std::move(coordinate));

        return std::forward<cell_t>(cell);
    }

    /*!\brief Tracks the last cell of a row within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `best_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a semi-global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_last_row_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_last_row_cell && !test_every_cell)
            invoke_comparator(cell, std::move(coordinate));
    }
    /*!\brief Tracks the last cell of a column within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `best_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a semi-global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_last_column_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_last_column_cell && !test_every_cell)
            invoke_comparator(cell, std::move(coordinate));

        return std::forward<cell_t>(cell);
    }

    /*!\brief Tracks the final cell of the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `best_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_final_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (!(test_every_cell || test_last_row_cell || test_last_column_cell))
            invoke_comparator(cell, std::move(coordinate));
    }
    //!\brief Resets the optimum such that a new alignment can be computed.
    void reset_optimum() noexcept
    {
        optimal_score = std::numeric_limits<score_type>::lowest();
        optimal_coordinate = {};
    }

    /*!\brief Handles the invocation of the optimum comparator and updater.
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `best_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \details
     *
     * Forwards the score and coordinate pair as a tuple and invokes the compare and set operation with the
     * so far best score/coordinate pair and the current score/coordinate pair.
     */
    template <typename cell_t>
    void invoke_comparator(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        compare_and_set_optimum(optimal_score, optimal_coordinate, cell.best_score(), std::move(coordinate));
    }

    /*!\brief Initialises the tracker with the row and column dimensions of the alignment matrix.
     *
     * \param[in] row_dimension The dimension of the rows, i.e. the number of rows.
     * \param[in] column_dimension The dimension of the columns, i.e. the number of columns.
     *
     * \details
     *
     * If the tracker uses an score update function specific for banded alignments, then the respective banded score
     * updater is initialised with the given matrix dimensions. Otherwise, this function is a noop and does nothing.
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Depends on the initialisation function of the used score updater function.
     */
    template <std::integral row_index_t, std::integral column_index_t>
    void initialise_tracker([[maybe_unused]] row_index_type<row_index_t> const row_dimension,
                            [[maybe_unused]] column_index_type<column_index_t> const column_dimension) noexcept
    {
        if constexpr (std::same_as<optimum_updater_t, max_score_banded_updater>)
            compare_and_set_optimum.set_target_indices(row_dimension, column_dimension);
    }
};
} // namespace seqan3::detail
