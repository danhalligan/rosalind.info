# Counting Optimal Alignments

from .helpers import Parser


def ctea(s1, s2):
    """Counting Optimal Alignments"""

    score, routes = {}, {}
    for j in range(len(s2) + 1):
        score[j, 0] = j
        routes[j, 0] = 1

    for i in range(len(s1) + 1):
        score[0, i] = i
        routes[0, i] = 1

    for j in range(len(s2)):
        for i in range(len(s1)):
            pos = [(j + 1, i), (j, i), (j, i + 1)]
            cost = [1, int(s1[i] != s2[j]), 1]
            scores = [score[pos[x]] + cost[x] for x in range(3)]
            best = min(scores)
            new = (j + 1, i + 1)
            score[new] = best
            routes[new] = sum(routes[pos[x]] for x in range(3) if scores[x] == best)

    return routes[len(s2), len(s1)] % 134217727


def main(file):
    seqs = Parser(file).seqs()
    print(ctea(seqs[0], seqs[1]))
