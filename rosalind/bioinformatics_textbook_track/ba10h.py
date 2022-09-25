# Estimate the Parameters of an HMM

from collections import defaultdict
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


def as_dict(x, r, c):
    g = defaultdict(float)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            g[r[i], c[j]] = x[i][j]
    return g


def estimate_tmat(path, states, to_dict=False):
    tmat = np.zeros((len(states), len(states)), dtype=float)
    for a, b in zip(path, path[1:]):
        tmat[states.index(a)][states.index(b)] += 1
    tmat = normalise(tmat, inc_zeros=True, min_val=1e-16)
    if to_dict:
        return as_dict(tmat, states, states)
    else:
        return tmat


def estimate_emat(seq, alphabet, path, states, to_dict=False):
    emat = np.zeros((len(states), len(alphabet)), dtype=float)
    for a, b in zip(path, seq):
        emat[states.index(a)][alphabet.index(b)] += 1
    emat = normalise(emat, inc_zeros=True, min_val=1e-16)
    if to_dict:
        return as_dict(emat, states, alphabet)
    else:
        return emat


def main(file):
    seq, alphabet, path, states = parse_input(open(file))
    tm = estimate_tmat(path, states)
    em = estimate_emat(seq, alphabet, path, states)
    print_mat(tm, states, states)
    print("--------")
    print_mat(em, states, alphabet)
