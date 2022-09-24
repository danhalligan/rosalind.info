# Estimate the Parameters of an HMM

import numpy as np
from .ba10e import normalise, print_mat


def parse_input(handle):
    seq = next(handle).rstrip()
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    path = next(handle).rstrip()
    next(handle)
    states = next(handle).split()
    return seq, alphabet, path, states


def tmat(path, states):
    tmat = np.zeros((len(states), len(states)), dtype=float)
    for a, b in zip(path, path[1:]):
        tmat[states.index(a)][states.index(b)] += 1
    return normalise(tmat, inc_zeros=True)


def emat(seq, alphabet, path, states):
    emat = np.zeros((len(states), len(alphabet)), dtype=float)
    for a, b in zip(path, seq):
        emat[states.index(a)][alphabet.index(b)] += 1
    return normalise(emat, inc_zeros=True)


def main(file):
    seq, alphabet, path, states = parse_input(open(file))
    tm = tmat(path, states)
    em = emat(seq, alphabet, path, states)
    print_mat(tm, states, states)
    print("--------")
    print_mat(em, states, alphabet)
