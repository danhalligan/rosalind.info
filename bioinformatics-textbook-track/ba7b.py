# Compute Limb Lengths in a Tree

import sys


def parse_mat(lines):
    return [[int(x) for x in y.split()] for y in lines]


def limb(j, x):
    limb_len = 100000
    n = len(x)
    for k in range(int(n)):
        for i in range(int(n)):
            if j != k and i != j:
                limb_len = min((x[i][j] + x[j][k] - x[i][k]) // 2, limb_len)
    return limb_len


def main(file):
    _, j, *arr = open(file).read().splitlines()
    print(limb(int(j), parse_mat(arr)))


if __name__ == "__main__":
    main(sys.argv[1])
