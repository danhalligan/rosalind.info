# Implement FarthestFirstTraversal

from math import sqrt


def read_types(handle, type):
    for line in handle:
        yield list(map(type, line.split()))


def euclidean_distance(a, b):
    """Euclidean distance between a pair of n-dimensional points"""
    return sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def ncd(x, centers):
    """Euclidean distance from DataPoint to its closest center"""
    return min(euclidean_distance(x, c) for c in centers)


def farthest_first_traversal(points, k):
    centers = [points[0]]
    while len(centers) < k:
        dists = [(i, ncd(point, centers)) for i, point in enumerate(points)]
        centers += [points[max(dists, key=lambda x: x[1])[0]]]
    return centers


def main(file):
    handle = open(file)
    k, m = next(read_types(handle, int))
    points = [point for point in read_types(handle, float)]
    centers = farthest_first_traversal(points, k)
    for c in centers:
        print(*c)
