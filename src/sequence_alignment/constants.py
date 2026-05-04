"""Shared scoring constants for DNA sequence alignment."""

GAP_PENALTY = 30

MISMATCH_COSTS = {
    ("A", "A"): 0,
    ("A", "C"): 110,
    ("A", "G"): 48,
    ("A", "T"): 94,
    ("C", "A"): 110,
    ("C", "C"): 0,
    ("C", "G"): 118,
    ("C", "T"): 48,
    ("G", "A"): 48,
    ("G", "C"): 118,
    ("G", "G"): 0,
    ("G", "T"): 110,
    ("T", "A"): 94,
    ("T", "C"): 48,
    ("T", "G"): 110,
    ("T", "T"): 0,
}

VALID_BASES = frozenset("ACGT")
