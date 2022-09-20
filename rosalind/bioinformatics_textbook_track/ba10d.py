# Compute the Probability of a String Emitted by an HMM

import numpy as np
from .ba10c import parse_input


def likelihood(seq, states, tmat, emat):
    mat = np.ones((len(seq) + 1, len(states)))

    for i, state in enumerate(states):
        mat[0, i] = emat[state, seq[0]] / len(states)

    for i, emission in enumerate(seq[1:], start=1):
        for j, state in enumerate(states):
            mat[i, j] = sum(
                tmat[prev, state] * emat[state, emission] * mat[i - 1, k]
                for k, prev in enumerate(states)
            )

    return sum(mat[i, :])


def main(file):
    seq, states, tmat, emat = parse_input(open(file))
    print(likelihood(seq, states, tmat, emat))
