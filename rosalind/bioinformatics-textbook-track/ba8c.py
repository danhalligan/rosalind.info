# Implement the Lloyd Algorithm for k-Means Clustering

from .ba8a import read_types, euclidean_distance
from functools import reduce


def ncd_assignment(x, centers):
    """Center index that minimises Euclidean distance to point"""
    dists = [euclidean_distance(x, c) for c in centers]
    return dists.index(min(dists))


def average_points(points):
    points = list(points)
    tot = reduce(lambda a, b: [a + b for a, b in zip(a, b)], points)
    return [x / len(points) for x in tot]


def compute_center(points, assigns, i):
    return average_points(p for p, a in zip(points, assigns) if a == i)


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
    points = [point for point in read_types(handle, float)]
    for m in k_means(points, k):
        print(*[f"{x:.3f}" for x in m])
