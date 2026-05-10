"""Memory-efficient sequence alignment implementation."""

import sys
import os
import time
import psutil
import threading
from dataclasses import dataclass
from pathlib import Path

# constants for DNA sequence alignment.
GAP_PENALTY = 30
MISMATCH_COSTS = {
    ("A", "A"): 0, ("A", "C"): 110, ("A", "G"): 48,
    ("A", "T"): 94, ("C", "A"): 110, ("C", "C"): 0,
    ("C", "G"): 118, ("C", "T"): 48, ("G", "A"): 48,
    ("G", "C"): 118, ("G", "G"): 0, ("G", "T"): 110,
    ("T", "A"): 94, ("T", "C"): 48, ("T", "G"): 110,
    ("T", "T"): 0,
}
VALID_BASES = frozenset("ACGT")

# data model
@dataclass(frozen=True)
class AlignmentResult:
    cost: int
    aligned_x: str
    aligned_y: str

# Input Parser
class InputFormatError(ValueError):
    """Raised when an input file does not match the expected project format."""

def _validate_base_sequence(sequence, label):
    if not sequence:
        raise InputFormatError(f"{label} sequence cannot be empty.")
    invalid_bases = sorted(set(sequence) - VALID_BASES)
    if invalid_bases:
        invalid_text = "".join(invalid_bases)
        raise InputFormatError(
            f"{label} sequence contains invalid DNA character(s): {invalid_text!r}."
        )
    return sequence

def _looks_like_integer(token):
    if token.isdecimal():
        return True
    return len(token) > 1 and token[0] in "+-" and token[1:].isdecimal()

def _parse_index_token(token):
    if not token.isdecimal():
        raise InputFormatError(
            f"Invalid insertion index {token!r}; indices must be non-negative integers."
        )
    return int(token)

def _validate_insertion_index(index, sequence_length, position):
    if isinstance(index, bool) or not isinstance(index, int):
        raise InputFormatError(
            f"Insertion index #{position + 1} must be an integer, got {index!r}."
        )
    if index < 0 or index >= sequence_length:
        raise InputFormatError(
            f"Insertion index #{position + 1} ({index}) is out of range for "
            f"sequence length {sequence_length}."
        )

def generate_sequence(base_sequence, insertion_indices):
    """Expand a base sequence using the project's insertion rules."""
    sequence = _validate_base_sequence(base_sequence, "Base")
    for position, insertion_index in enumerate(insertion_indices):
        _validate_insertion_index(insertion_index, len(sequence), position)
        sequence = (
            sequence[: insertion_index + 1]
            + sequence
            + sequence[insertion_index + 1 :]
        )
    return sequence

def parse_problem_lines(lines):
    """Parse input lines and return the generated sequence pair."""
    tokens = [line.strip() for line in lines if line.strip()]
    if not tokens:
        raise InputFormatError("Input file is empty.")

    if _looks_like_integer(tokens[0]):
        raise InputFormatError("Expected first base sequence before insertion indices.")

    first_base = _validate_base_sequence(tokens[0], "First base")
    first_indices = []
    second_base = None
    second_indices = []

    for token in tokens[1:]:
        if _looks_like_integer(token):
            index = _parse_index_token(token)
            if second_base is None:
                first_indices.append(index)
            else:
                second_indices.append(index)
            continue

        if second_base is None:
            second_base = _validate_base_sequence(token, "Second base")
            continue

        raise InputFormatError(f"Unexpected token after second sequence: {token!r}.")

    if second_base is None:
        raise InputFormatError("Expected second base sequence after first insertion indices.")

    return (
        generate_sequence(first_base, first_indices),
        generate_sequence(second_base, second_indices),
    )

def read_problem(input_path):
    path = Path(input_path)
    with path.open("r", encoding="utf-8") as input_file:
        return parse_problem_lines(input_file)


# Metrics
@dataclass(frozen=True)
class RunMetrics:
    elapsed_seconds: float
    memory_before_bytes: int
    memory_after_bytes: int
    peak_memory_bytes: int

    @property
    def elapsed_ms(self):
        return self.elapsed_seconds * 1000

    @property
    def peak_memory_kb(self):
        return self.peak_memory_bytes / 1024

def run_with_metrics(operation, sample_interval_seconds=0.001):
    process = psutil.Process(os.getpid())
    memory_before = process.memory_info().rss
    peak_memory = memory_before
    stop_event = threading.Event()

    def sample_memory():
        nonlocal peak_memory
        while not stop_event.is_set():
            peak_memory = max(peak_memory, process.memory_info().rss)
            stop_event.wait(sample_interval_seconds)

    sampler = threading.Thread(target=sample_memory, daemon=True)
    start_time = time.perf_counter()
    sampler.start()
    try:
        value = operation()
    finally:
        elapsed_seconds = time.perf_counter() - start_time
        stop_event.set()
        sampler.join()

    memory_after = process.memory_info().rss
    peak_memory = max(peak_memory, memory_after)
    return value, RunMetrics(elapsed_seconds, memory_before, memory_after, peak_memory)


# Algorithm
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


# Main (Connect everything)
def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    x,y = read_problem(input_file)

    result, metrics = run_with_metrics(lambda:align_efficient(x,y))

    with open(output_file, "w") as f:
        f.write(str(result.cost) + "\n")
        f.write(result.aligned_x + "\n")
        f.write(result.aligned_y + "\n")
        f.write(f"{metrics.elapsed_ms}\n")
        f.write(f"{metrics.peak_memory_kb}\n")

if __name__ == "__main__":
    main()