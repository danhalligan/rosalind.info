# Implement the Neighbor Joining Algorithm

from collections import defaultdict
import numpy as np
from .ba7b import parse_mat
from .ba7d import as_edges, closest


def nj_matrix(D, n):
    ND = np.copy(D)
    for i in range(len(D)):
        for j in range(len(D)):
            if i != j:
                ND[i, j] = (n - 2) * D[i, j] - sum(D[i, :]) - sum(D[j, :])
    return ND


def neighbor_joining(D, n, labels=None):
    if not labels:
        labels = list(range(n))

    if n == 2:
        T = defaultdict(list)
        T[labels[0]].append({"n": labels[1], "w": D[0][1]})
        return T

    ND = nj_matrix(D, n)
    i, j = closest(ND)
    delta = (sum(D[i, :]) - sum(D[j, :])) / (n - 2)
    limb_i = (D[i, j] + delta) / 2
    limb_j = (D[i, j] - delta) / 2

    li = labels[i]
    lj = labels[j]

    D = np.append(D, np.zeros((1, len(D))), axis=0)
    D = np.append(D, np.zeros((len(D), 1)), axis=1)
    labels = labels + [max(labels) + 1]

    for k in range(n):
        D[k, n] = (D[k, i] + D[k, j] - D[i, j]) / 2
        D[n, k] = (D[k, i] + D[k, j] - D[i, j]) / 2
    for x in [j, i]:
        D = np.delete(D, x, 0)
        D = np.delete(D, x, 1)
        del labels[x]

    T = neighbor_joining(D, n - 1, labels)

    T[labels[-1]].append({"n": li, "w": limb_i})
    T[labels[-1]].append({"n": lj, "w": limb_j})
    return T


def main(file):
    n, *D = open(file).read().splitlines()
    D = np.array(parse_mat(D), float)
    graph = neighbor_joining(D, int(n))
    for edge in as_edges(graph):
        print(edge)
