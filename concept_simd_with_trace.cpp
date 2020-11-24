

// traceback matrix with only one value?

// inside the recursion we need to decide for the best trace
// we need to know the best trace to the next cell.
// right now we propagate the cell value to the new function

// new score = score + diagonal
// we already know from which direction the trace came
// we have to always load an entire column even though they are not related.
// can I turn around the computation?

// Score matrix: value_type
// Also this is a proxy type that we return.

// In case of the trace matrix we need another score matrix.
// pair<score_t, score_t> -> best_left + left
// pair<score_t, score_t> -> best_up + up  // One additional value in the vertical column

// So it is aware of the difference.
// previous_cell layout: <best_left, left>, <best_up, up>

// compute_inner_cell(diagonal, cell, this->scoring_scheme.score(alphabet1, alphabet2))
template <typename score_t, typename previous_entries_t, typename score_t>
auto compute_cell(score_t best_score,
                  previous_entries_t previous,
                  score_t substitution_score)
{
    auto && [previous_left, previous_up] = previous;

    // Diagonal value + substitution score.
    best_score += substitution_score;

    // Compute the optimum coming from the upper cell
    previous_up.best_score() += gap_open;
    previous_up.up_score() += gap_extension;

    auto cmp_mask = (previous_up.best_score() < previous_up.up_score());
    previous_up.up_score() = cmp_mask ? previous_up.up_score() : previous_up.best_score();
    trace_t trace_up = cmp_mask ? trace_directions::up : trace_directions::up_open;

    // Maximum between U and D.
    cmp_mask = (best_score < previous_up.up_score())
    best_score = cmp_mask ? previous_up.up_score() : best_score;
    trace_t trace_best = cmp_mask ? trace_up : trace_directions::diagonal | trace_up;

    // Compute the optimum coming from the left cell
    substitution_score = previous_left.best_score() + gap_open;
    previous_up.best_score() = previous_left.left_score() + gap_extension;

    cmp_mask = (previous_up.best_score() < substitution_score);
    substitution_score = cmp_mask ? substitution_score : previous_up.best_score();
    trace_up |= cmp_mask ? trace_directions::left : trace_directions::left_open;

    // Maximum between L and D.
    // Order of this has influence of precedence in trace back.
    // So here it is better to come from left than from up or diagonal.
    // This means if we get into this cell, we need to follow left instead
    // of up.
    // In a clean interface this would be coupled.
    // So the same structure that decides how to handle this must also
    // decide how to follow the trace. Or we need to document it carefully.
    // So the aligned sequence builder and this operation.
    cmp_mask = (best_score < substitution_score)
    best_score = cmp_mask ? substitution_score : best_score;
    trace_best = cmp_mask ? trace_up : trace_best | trace_left;

    // Now we have everything together!
    return {{best_score, score}, {best_score, previous_up.up_score()}, trace_best};

}

compute_cell(diagonal, previous_cell)
score = last_score + diagonal;

tmp = previous_left.best_score() + gap_open; // We need this already.
score_from_left = previous_cell.horizontal_score() + gap_extension;
score_from_left = (score_from_left < tmp) ? tmp : score_from_left;

// check trace value here

tmp = previous_up.best_score() + gap_open;
score_from_above = previous_cell.up_score() + gap_extension;
score_from_above = (score_from_above < tmp) ? tmp : score_from_above;

// check trace value here

score = max(score_from_left, score_from_above, score);
return {score, score_from_left, score_from_above};
// we would only need one trace value to store!

