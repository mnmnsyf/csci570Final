"""Shared result model for sequence alignment implementations."""

from dataclasses import dataclass


@dataclass(frozen=True)
class AlignmentResult:
    cost: int
    aligned_x: str
    aligned_y: str
