# Compute the Probability of a Hidden Path

from math import prod


def ba10a(seq, states, tmat):
    return 0.5 * prod(tmat[x] for x in zip(seq, seq[1:]))


def parse_input(handle):
    seq = next(handle).rstrip()
    next(handle)
    states = next(handle).split()
    next(handle)
    tmat = dict(
        ((states[i], states[j]), float(v))
        for i, x in enumerate(handle.read().splitlines()[1:])
        for j, v in enumerate(x.split()[1:])
    )
    return seq, states, tmat


def main(file):
    pr = ba10a(*parse_input(open(file)))
    print(f"{pr:.11e}")
