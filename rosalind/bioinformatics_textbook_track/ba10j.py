# Solve the Soft Decoding Problem

import numpy as np
from .ba10c import parse_input


def forward(seq, states, tmat, emat):
    mat = np.ones((len(seq), len(states)))

    for i, state in enumerate(states):
        mat[0, i] = emat[state, seq[0]]
    for i, emission in enumerate(seq[1:], start=1):
        for j, state in enumerate(states):
            mat[i, j] = sum(
                tmat[prev, state] * emat[state, emission] * mat[i - 1, k]
                for k, prev in enumerate(states)
            )

    return mat


def backward(seq, states, tmat, emat):
    mat = np.ones((len(seq), len(states)))

    for i, emission in enumerate(seq[::-1][:-1], start=1):
        for j, state in enumerate(states):
            mat[len(seq) - i - 1, j] = sum(
                tmat[state, prev] * emat[prev, emission] * mat[len(seq) - i, k]
                for k, prev in enumerate(states)
            )
    return mat


def soft_decode(seq, states, tmat, emat, normalise=True):
    tot = forward(seq, states, tmat, emat) * backward(seq, states, tmat, emat)
    if normalise:
        tot = tot / np.sum(tot, axis=1, keepdims=True)
    return tot


def main(file):
    seq, states, tmat, emat = parse_input(open(file))
    tot = soft_decode(seq, states, tmat, emat)
    print(*states, sep="\t")
    for r in np.round(tot, 4):
        print(*r, sep="\t")
