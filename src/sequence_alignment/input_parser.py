"""Skeleton parser for the CSCI570 sequence-generation input format."""

from pathlib import Path


class InputFormatError(ValueError):
    """Raised when an input file does not match the expected project format."""


def generate_sequence(base_sequence, insertion_indices):
    """TODO: expand a base sequence using the project's insertion rules."""
    return base_sequence


def parse_problem_lines(lines):
    """TODO: parse input lines and return the generated sequence pair."""
    tokens = [line.strip() for line in lines if line.strip()]
    if not tokens:
        raise InputFormatError("Input file is empty.")
    return "", ""


def read_problem(input_path):
    path = Path(input_path)
    with path.open("r", encoding="utf-8") as input_file:
        return parse_problem_lines(input_file)
