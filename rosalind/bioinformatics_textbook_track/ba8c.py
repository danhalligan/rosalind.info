# Implement the Lloyd Algorithm for k-Means Clustering

import numpy as np
from .ba8a import read_types, euclidean_distance


def ncd_assignment(x, centers):
    """Center index that minimises Euclidean distance to point"""
    dists = [euclidean_distance(x, c) for c in centers]
    return dists.index(min(dists))


def compute_center(points, assigns, i):
    data = [p for p, a in zip(points, assigns) if a == i]
    return np.mean(np.array(data), 0)


# I'm lazily avoiding checking for convergence here and hoping that 20
# iterations is enough. This converges quite quickly though...
def k_means(points, k, iter=20):
    centers = points[:k]
    for _ in range(iter):
        assigns = [ncd_assignment(point, centers) for point in points]
        centers = [compute_center(points, assigns, i) for i in range(k)]
    return centers


def main(file):
    handle = open(file)
    k, m = next(read_types(handle, int))
    points = [np.array(point) for point in read_types(handle, float)]
    for m in k_means(points, k):
        print(*[f"{x:.3f}" for x in m])
