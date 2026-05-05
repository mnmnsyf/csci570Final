"""Memory-efficient sequence alignment implementation."""

from sequence_alignment.alignment import AlignmentResult
from sequence_alignment.constants import GAP_PENALTY, MISMATCH_COSTS


def _character_cost(char_x, char_y):
    return MISMATCH_COSTS[(char_x, char_y)]


def _alignment_cost(aligned_x, aligned_y):
    total_cost = 0
    for char_x, char_y in zip(aligned_x, aligned_y):
        if char_x == "_" or char_y == "_":
            total_cost += GAP_PENALTY
        else:
            total_cost += _character_cost(char_x, char_y)
    return total_cost


def linear_space_scores(sequence_x, sequence_y):
    """Return the final DP row using O(len(sequence_y)) memory."""
    # Base row: aligning an empty first sequence to prefixes of sequence_y.
    previous_row = [column * GAP_PENALTY for column in range(len(sequence_y) + 1)]

    for row, char_x in enumerate(sequence_x, start=1):
        # Only the previous row is needed to build the current row.
        current_row = [row * GAP_PENALTY]
        for column, char_y in enumerate(sequence_y, start=1):
            diagonal_cost = previous_row[column - 1] + _character_cost(char_x, char_y)
            delete_cost = previous_row[column] + GAP_PENALTY
            insert_cost = current_row[column - 1] + GAP_PENALTY
            current_row.append(min(diagonal_cost, delete_cost, insert_cost))
        previous_row = current_row

    return previous_row


def _align_with_table(sequence_x, sequence_y):
    # Hirschberg uses this full-table DP only for tiny base cases.
    rows = len(sequence_x) + 1
    columns = len(sequence_y) + 1
    cost_table = [[0] * columns for _ in range(rows)]

    for row in range(1, rows):
        cost_table[row][0] = row * GAP_PENALTY
    for column in range(1, columns):
        cost_table[0][column] = column * GAP_PENALTY

    for row, char_x in enumerate(sequence_x, start=1):
        for column, char_y in enumerate(sequence_y, start=1):
            diagonal_cost = cost_table[row - 1][column - 1] + _character_cost(char_x, char_y)
            delete_cost = cost_table[row - 1][column] + GAP_PENALTY
            insert_cost = cost_table[row][column - 1] + GAP_PENALTY
            cost_table[row][column] = min(diagonal_cost, delete_cost, insert_cost)

    aligned_x = []
    aligned_y = []
    row = len(sequence_x)
    column = len(sequence_y)

    while row > 0 or column > 0:
        if row > 0 and column > 0:
            diagonal_cost = cost_table[row - 1][column - 1] + _character_cost(
                sequence_x[row - 1], sequence_y[column - 1]
            )
            if cost_table[row][column] == diagonal_cost:
                aligned_x.append(sequence_x[row - 1])
                aligned_y.append(sequence_y[column - 1])
                row -= 1
                column -= 1
                continue

        if row > 0 and cost_table[row][column] == cost_table[row - 1][column] + GAP_PENALTY:
            aligned_x.append(sequence_x[row - 1])
            aligned_y.append("_")
            row -= 1
        else:
            aligned_x.append("_")
            aligned_y.append(sequence_y[column - 1])
            column -= 1

    return "".join(reversed(aligned_x)), "".join(reversed(aligned_y))


def _find_split(sequence_x, sequence_y, middle):
    left_scores = linear_space_scores(sequence_x[:middle], sequence_y)

    # Reverse the suffixes so each score lines up with a possible split in sequence_y.
    right_scores = linear_space_scores(sequence_x[middle:][::-1], sequence_y[::-1])

    best_column = 0
    best_cost = left_scores[0] + right_scores[-1]
    last_column = len(sequence_y)

    for column in range(1, last_column + 1):
        split_cost = left_scores[column] + right_scores[last_column - column]
        if split_cost < best_cost:
            best_cost = split_cost
            best_column = column

    return best_column


def _hirschberg(sequence_x, sequence_y):
    if not sequence_x:
        return "_" * len(sequence_y), sequence_y
    if not sequence_y:
        return sequence_x, "_" * len(sequence_x)
    if len(sequence_x) == 1 or len(sequence_y) == 1:
        return _align_with_table(sequence_x, sequence_y)

    # Split sequence_x in half and choose the sequence_y split with minimum combined cost.
    middle = len(sequence_x) // 2
    split_column = _find_split(sequence_x, sequence_y, middle)

    left_x, left_y = _hirschberg(sequence_x[:middle], sequence_y[:split_column])
    right_x, right_y = _hirschberg(sequence_x[middle:], sequence_y[split_column:])
    return left_x + right_x, left_y + right_y


def align_efficient(sequence_x, sequence_y):
    """Align two sequences using Hirschberg's divide-and-conquer algorithm."""
    aligned_x, aligned_y = _hirschberg(sequence_x, sequence_y)
    return AlignmentResult(_alignment_cost(aligned_x, aligned_y), aligned_x, aligned_y)
