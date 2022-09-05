# Compute the Number of Breakpoints in a Permutation

from .ba6a import read_perm


def nbreakpoints(perm):
    perm = [0] + perm + [len(perm) + 1]
    return sum(perm[i] + 1 != perm[i + 1] for i in range(len(perm) - 1))


def main(file):
    s = open(file).read().rstrip()
    print(nbreakpoints(read_perm(s)))
