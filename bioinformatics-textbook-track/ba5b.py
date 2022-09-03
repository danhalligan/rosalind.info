# Find the Length of a Longest Path in a Manhattan-like Grid

import numpy as np


def manhattan_tourist(n, m, down, right):
    s = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max([s[i - 1][j] + down[i - 1][j], s[i][j - 1] + right[i][j - 1]])
    return int(s[n][m])


def main(file):
    dim, *lines = open(file).read().splitlines()
    n, m = map(int, dim.split())
    down = np.array([list(map(int, x.split())) for x in lines[0:n]])
    right = np.array([list(map(int, x.split())) for x in lines[(n + 1) :]])
    print(manhattan_tourist(n, m, down, right))
