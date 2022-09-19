# Implement the Viterbi Algorithm

from math import log
import numpy as np


def parse_input(handle):
    seq = next(handle).rstrip()
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    states = next(handle).split()
    next(handle)
    lines = [next(handle) for _ in range(len(states) + 1)]
    tmat = {
        (states[i], states[j]): log(float(v))
        for i, x in enumerate(lines[1:])
        for j, v in enumerate(x.split()[1:])
    }
    next(handle)
    lines = [next(handle) for i in range(len(states) + 1)]
    emat = {
        (states[i], alphabet[j]): log(float(v))
        for i, x in enumerate(lines[1:])
        for j, v in enumerate(x.split()[1:])
    }
    return seq, states, tmat, emat


def viterbi(seq, states, tmat, emat):
    mat = np.zeros((len(seq) + 1, len(states)))
    ptr = np.zeros((len(seq) + 1, len(states)), dtype=int)
    for i, emission in enumerate(seq, start=1):
        for j, state in enumerate(states):
            opt = [
                tmat[prev, state] + emat[state, emission] + mat[i - 1, k]
                for k, prev in enumerate(states)
            ]
            p = opt.index(max(opt))
            ptr[i, j] = p
            mat[i, j] = max(opt)
    ind = np.argmax(mat[i, :])
    seq = states[ind]
    while i > 1:
        seq = states[ptr[i, ind]] + seq
        ind = ptr[i, ind]
        i -= 1
    return seq


def main(file):
    seq, states, tmat, emat = parse_input(open(file))
    print(viterbi(seq, states, tmat, emat))
