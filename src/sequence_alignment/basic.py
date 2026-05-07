"""Basic O(mn)-time and O(mn)-space DP sequence alignment."""

from sequence_alignment.alignment import AlignmentResult
from sequence_alignment.constants import GAP_PENALTY, MISMATCH_COSTS


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
