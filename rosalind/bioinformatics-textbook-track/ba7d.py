# Implement UPGMA

from .ba7b import parse_mat
import numpy as np


def as_edges(graph):
    edges = []
    for k in sorted(graph):
        for v in graph[k]:
            edges += [f"{k}->{v['n']}:{v['w']:.3f}"]
            edges += [f"{v['n']}->{k}:{v['w']:.3f}"]
    return sorted(edges)


# find (first) minimum off diagonal index in an array
def closest(D):
    D = np.copy(D)
    np.fill_diagonal(D, D.max() + 1)
    return divmod(D.argmin(), D.shape[1])


# replace the ith row/col with the average of the ith and jth and remove
# the jth
def average_ind(D, i, j, di, dj):
    D = np.copy(D)
    av = (
        D[
            i,
        ]
        * di
        + D[
            j,
        ]
        * dj
    ) / (di + dj)
    D[i, :] = av
    D[:, i] = av
    D = np.delete(D, j, 0)
    D = np.delete(D, j, 1)
    np.fill_diagonal(D, 0)
    return D


def upgma(D, m):
    clusters = list(range(0, m))
    ages = {}
    deg = {}
    for x in clusters:
        ages[x] = 0
        deg[x] = 1
    T = {}

    # a label for internal nodes as we add them
    node = m
    while len(clusters) > 1:
        # find the two closest clusters Ci and Cj
        i, j = closest(D)
        a, b = clusters[i], clusters[j]

        T[node] = [
            {"n": a, "w": D[i, j] / 2 - ages[a]},
            {"n": b, "w": D[i, j] / 2 - ages[b]},
        ]
        deg[node] = deg[a] + deg[b]
        ages[node] = D[i, j] / 2
        clusters[i] = node
        del clusters[j]
        D = average_ind(D, *closest(D), deg[a], deg[b])
        node += 1

    return T


def main(file):
    n, *D = open(file).read().splitlines()
    D = np.array(parse_mat(D), float)
    graph = upgma(D, int(n))
    for edge in as_edges(graph):
        print(edge)
