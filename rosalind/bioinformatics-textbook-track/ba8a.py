# Implement FarthestFirstTraversal

from math import sqrt


def euclidean_distance(a, b):
    """Euclidean distance between a pair of n-dimensional points"""
    return sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def nearest_center_distance(point, centers):
    """Euclidean distance from DataPoint to its closest center"""
    return min(euclidean_distance(point, center) for center in centers)


def farthest_first_traversal(points, k):
    centers = [points[0]]
    while len(centers) < k:
        dists = [
            (i, nearest_center_distance(point, centers))
            for i, point in enumerate(points)
        ]
        centers += [points[max(dists, key=lambda x: x[1])[0]]]
    return centers


def main(file):
    info, *points = open(file).read().splitlines()
    k, m = map(int, info.split())
    points = [list(map(float, x.split())) for x in points]
    centers = farthest_first_traversal(points, k)
    for c in centers:
        print(*c)
