# Local Alignment with Scoring Matrix and Affine Gap Penalty

from .helpers import Parser
from rosalind.helpers import blosum62


def laff(v, w, score, sigma, epsilon):
    """Local Alignment with Affine Gap Penalty"""

    # Initialize
    x = [[0 for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
    m = [[0 for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
    y = [[0 for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
    backtrack = [[0 for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
    max_score = -1
    max_i, max_j = 0, 0

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            x[i][j] = max([x[i - 1][j] + epsilon, m[i - 1][j] + sigma])
            y[i][j] = max([y[i][j - 1] + epsilon, m[i][j - 1] + sigma])
            s = [x[i][j], m[i - 1][j - 1] + score[v[i - 1]][w[j - 1]], y[i][j], 0]
            m[i][j] = max(s)
            backtrack[i][j] = s.index(m[i][j])

            if m[i][j] > max_score:
                max_score = m[i][j]
                max_i, max_j = i, j

    # Backtrack to start of the local alignment
    # Recover substrings (alignment not required)
    i, j = max_i, max_j
    va, wa = v[:i], w[:j]
    while backtrack[i][j] != 3 and i * j != 0:
        if backtrack[i][j] == 0:
            i -= 1
        elif backtrack[i][j] == 1:
            i -= 1
            j -= 1
        elif backtrack[i][j] == 2:
            j -= 1
    va, wa = va[i:], wa[j:]

    return {"score": str(max_score), "va": va, "wa": wa}


def main(file):
    seqs = Parser(file).seqs()
    res = laff(seqs[0], seqs[1], blosum62(), -11, -1)
    print(res["score"])
    print(res["va"])
    print(res["wa"])
