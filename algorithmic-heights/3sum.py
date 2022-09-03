# 3SUM

# Naively we would check all subsequent elements for each element
# and terminate on finding a valid pair (O[n^2]).
# We'll use a hash for O(n)

from helpers import ints


def three_sum(n, a):
    h = {}
    for i in range(n):
        for j in range(i + 1, n):
            s = a[i] + a[j]
            if s in h:
                return [h[s] + 1, i + 1, j + 1]
        h[-a[i]] = i
    return [-1]


def main(file):
    info, *arrays = open(file).read().splitlines()
    k, n = ints(info)
    for arr in arrays:
        print(*three_sum(int(n), ints(arr)))
