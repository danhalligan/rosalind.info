# Compute the Probability of an Outcome Given a Hidden Path

from math import prod


def parse_input(handle):
    seq = next(handle).rstrip()
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    path = next(handle).rstrip()
    next(handle)
    states = next(handle).split()
    next(handle)
    tmat = {
        (states[i], alphabet[j]): float(v)
        for i, x in enumerate(handle.read().splitlines()[1:])
        for j, v in enumerate(x.split()[1:])
    }
    return seq, path, tmat


def ba10b(seq, path, tmat):
    return prod(tmat[x] for x in zip(path, seq))


def main(file):
    pr = ba10b(*parse_input(open(file)))
    print(f"{pr:.11e}")
