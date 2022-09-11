# Compute the Squared Error Distortion

from .ba8a import ncd, read_types


def distortion(points, centers):
    return (1 / len(points)) * sum(ncd(point, centers) ** 2 for point in points)


def main(file):
    handle = open(file)
    k, _ = next(read_types(handle, int))
    gen = read_types(handle, float)
    centers = [next(gen) for i in range(k)]
    _ = next(handle)
    points = [point for point in gen]
    print(round(distortion(points, centers), 3))
