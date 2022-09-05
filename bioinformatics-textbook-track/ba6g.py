# Implement ChromosomeToCycle

import sys
from .ba6a import read_perm, format_perm


def cycle2chromosome(cycle):
    nodes = []
    for j1, j2 in zip(cycle[::2], cycle[1::2]):
        if j1 < j2:
            nodes += [j2 // 2]
        else:
            nodes += [-j1 // 2]
    return nodes


def ba6f(s):
    return cycle2chromosome(read_perm(s))


def main(file):
    s = open(file).read().rstrip()
    print(format_perm(ba6f(s)))


if __name__ == "__main__":
    main(sys.argv[1])
