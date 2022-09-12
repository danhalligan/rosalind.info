# Implement Hierarchical Clustering

import numpy as np
from collections import defaultdict
from .ba7d import closest, average_ind


def descendants(T, node):
    q = [node]
    x = []
    while len(q):
        n = q.pop(0)
        if n in T:
            q += T[n]
        else:
            x += [n]
    return x


def hierarchical_clustering(D, n):
    clusters = list(range(1, n + 1))
    T = {}
    size = defaultdict(lambda: 1)  # the number of descendants of a node
    node = n
    while len(clusters) > 1:
        node += 1
        i, j = closest(D)
        a, b = clusters[i], clusters[j]
        T[node] = [a, b]
        size[node] = size[a] + size[b]
        D = average_ind(D, *closest(D), size[a], size[b])
        clusters[i] = node
        del clusters[j]
        yield descendants(T, a) + descendants(T, b)


def main(file):
    n, *m = open(file).read().splitlines()
    m = np.array([list(map(float, x.split())) for x in m])
    for step in hierarchical_clustering(m, int(n)):
        print(*step)
