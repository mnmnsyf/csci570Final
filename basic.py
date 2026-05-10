"""Basic O(mn)-time and O(mn)-space DP sequence alignment."""

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
def build_cost_table(sequence_x, sequence_y):
    """Build and return the full (m+1) x (n+1) DP cost table."""
    m, n = len(sequence_x), len(sequence_y)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        dp[i][0] = i * GAP_PENALTY
    for j in range(n + 1):
        dp[0][j] = j * GAP_PENALTY

    for i in range(1, m + 1):
        xi = sequence_x[i - 1]
        for j in range(1, n + 1):
            dp[i][j] = min(
                MISMATCH_COSTS[(xi, sequence_y[j - 1])] + dp[i - 1][j - 1],
                GAP_PENALTY + dp[i - 1][j],
                GAP_PENALTY + dp[i][j - 1],
            )
    return dp


def align_basic(sequence_x, sequence_y):
    """O(mn)-space sequence alignment. Returns AlignmentResult."""
    m, n = len(sequence_x), len(sequence_y)
    dp = build_cost_table(sequence_x, sequence_y)
    cost = dp[m][n]

    ax, ay = [], []
    i, j = m, n
    while i > 0 and j > 0:
        diag = MISMATCH_COSTS[(sequence_x[i - 1], sequence_y[j - 1])] + dp[i - 1][j - 1]
        if dp[i][j] == diag:
            ax.append(sequence_x[i - 1])
            ay.append(sequence_y[j - 1])
            i -= 1; j -= 1
        elif dp[i][j] == GAP_PENALTY + dp[i - 1][j]:
            ax.append(sequence_x[i - 1])
            ay.append('_')
            i -= 1
        else:
            ax.append('_')
            ay.append(sequence_y[j - 1])
            j -= 1

    while i > 0:
        ax.append(sequence_x[i - 1]); ay.append('_'); i -= 1
    while j > 0:
        ax.append('_'); ay.append(sequence_y[j - 1]); j -= 1

    return AlignmentResult(cost, ''.join(reversed(ax)), ''.join(reversed(ay)))

# Main (Connect everything)
def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    x,y = read_problem(input_file)

    result, metrics = run_with_metrics(lambda:align_basic(x,y))

    with open(output_file, "w") as f:
        f.write(str(result.cost) + "\n")
        f.write(result.aligned_x + "\n")
        f.write(result.aligned_y + "\n")
        f.write(f"{metrics.elapsed_ms}\n")
        f.write(f"{metrics.peak_memory_kb}\n")

if __name__ == "__main__":
    main()
