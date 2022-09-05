# 2-Way Partition

from .helpers import ints
import sys


def par(x):
    pivot, pi = x[0], 0
    for i in range(1, len(x)):
        if x[i] <= pivot:
            pi += 1
            x[i], x[pi] = x[pi], x[i]
    x[0], x[pi] = x[pi], x[0]
    return x


def main(file):
    n, x = open(file).read().splitlines()
    print(*par(ints(x)))


if __name__ == "__main__":
    main(sys.argv[1])
