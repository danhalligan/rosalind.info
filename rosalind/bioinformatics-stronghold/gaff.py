# Global Alignment with Scoring Matrix and Affine Gap Penalty

from .helpers import Parser
from rosalind.helpers import blosum62


def insert_indel(word, i):
    return word[:i] + "-" + word[i:]


def gaff(v, w, score, sigma, epsilon):
    """Global Alignment with Affine Gap Penalty"""
    # See Biological Sequence Analysis page 29
    # We now have to keep track of three matrices m, x and y

    m = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    x = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    y = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    px = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    pm = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    py = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]

    # Initialize the edges with the given penalties.
    for i in range(1, len(v) + 1):
        x[i][0] = sigma + (i - 1) * epsilon
        m[i][0] = sigma + (i - 1) * epsilon
        y[i][0] = 10 * sigma
    for j in range(1, len(w) + 1):
        y[0][j] = sigma + (j - 1) * epsilon
        m[0][j] = sigma + (j - 1) * epsilon
        x[0][j] = 10 * sigma

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s = [x[i - 1][j] + epsilon, m[i - 1][j] + sigma]
            x[i][j] = max(s)
            px[i][j] = s.index(x[i][j])

            s = [y[i][j - 1] + epsilon, m[i][j - 1] + sigma]
            y[i][j] = max(s)
            py[i][j] = s.index(y[i][j])

            s = [x[i][j], m[i - 1][j - 1] + score[v[i - 1]][w[j - 1]], y[i][j]]
            m[i][j] = max(s)
            pm[i][j] = s.index(m[i][j])

    # Initialize the values of i, j and the aligned sequences.
    i, j = len(v), len(w)
    a1, a2 = v, w

    # Get the maximum score, and the corresponding backtrack starting position.
    scores = [x[i][j], m[i][j], y[i][j]]
    max_score = max(scores)
    s = scores.index(max_score)

    # Backtrack to the edge of the matrix starting bottom right.
    while i * j != 0:
        if s == 0:
            if px[i][j] == 1:
                s = 1
            i -= 1
            a2 = insert_indel(a2, j)
        elif s == 1:
            if pm[i][j] == 1:
                i -= 1
                j -= 1
            else:
                s = pm[i][j]
        else:
            if py[i][j] == 1:
                s = 1
            j -= 1
            a1 = insert_indel(a1, i)

    # Prepend the necessary preceding indels to get to (0,0).
    for _ in range(i):
        a2 = insert_indel(a2, 0)
    for _ in range(j):
        a1 = insert_indel(a1, 0)

    return {"dist": str(max_score), "a1": a1, "a2": a2}


def main(file):
    seqs = Parser(file).seqs()
    res = gaff(seqs[0], seqs[1], blosum62(), -11, -1)
    print(res["dist"])
    print(res["a1"])
    print(res["a2"])
