# Find a Middle Edge in an Alignment Graph in Linear Space

from math import floor
from .ba5e import blosum62


# Calculate scores in linear space
# This keeps a vector of scores scores for the current and previous alignemnt
# columns (sc and psc). It also computes a single backtrack vector (bt)
# corresponding to cell that each value of sc was computed from.
def calculate_scores(s1, s2, scores, penalty):
    sc = list(range(0, (len(s1) + 1) * penalty, penalty))
    bt = [0] * (len(s1) + 1)
    for j in range(1, len(s2) + 1):
        psc = sc[:]
        sc[0] = psc[0] + penalty
        for i in range(1, len(s1) + 1):
            opt = [
                psc[i] + penalty,
                sc[i - 1] + penalty,
                psc[i - 1] + scores[s1[i - 1]][s2[j - 1]],
            ]
            sc[i] = max(opt)
            bt[i] = opt.index(sc[i])
    return sc, bt


# To find, the middle edge, we first find the highest scoring value of the
# middle column, then find the edge from this node to the next.
def middle_edge(s1, s2, scores, penalty):
    mid = floor((len(s2)) / 2)
    sc1, _ = calculate_scores(s1, s2[:mid], scores, penalty)
    sc2, bt2 = calculate_scores(s1[::-1], s2[mid::][::-1], scores, penalty)
    total = [a + b for (a, b) in zip(sc1, sc2[::-1])]
    best = total.index(max(total))
    n1 = (best, mid)
    moves = [(n1[0], n1[1] + 1), (n1[0] + 1, n1[1]), (n1[0] + 1, n1[1] + 1)]
    n2 = moves[bt2[::-1][best]]
    return (n1, n2)


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(middle_edge(s1, s2, blosum62(), -5))
