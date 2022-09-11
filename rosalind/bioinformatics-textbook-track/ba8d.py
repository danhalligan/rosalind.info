# Implement the Soft k-Means Clustering Algorithm

import numpy as np
from .ba8a import read_types, euclidean_distance


def partition_function(x, centers, β):
    num = [np.exp(-β * euclidean_distance(x, cen)) for cen in centers]
    return np.array(num) / sum(num)


def hidden_matrix(data, centers, β):
    return np.array([partition_function(x, centers, β) for x in data])


def soft_k_means(points, k, β, iter=20):
    centers = np.array(points[:k])
    points = np.array(points)
    for _ in range(iter):
        hm = hidden_matrix(points, centers, β)
        centers = [np.dot(hm[:, i], points) for i in range(k)]
        sums = np.sum(hm, 0)
        centers = np.transpose(np.transpose(centers) / sums)
    return centers


def main(file):
    handle = open(file)
    k, m = next(read_types(handle, int))
    β = next(read_types(handle, float))[0]
    points = [np.array(point) for point in read_types(handle, float)]
    for m in soft_k_means(points, k, β):
        print(*[f"{x:.3f}" for x in m])
