"""Parser for the CSCI570 sequence-generation input format."""

from pathlib import Path

from sequence_alignment.constants import VALID_BASES


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
