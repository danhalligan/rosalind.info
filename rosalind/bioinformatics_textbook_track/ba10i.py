# Implement Viterbi Learning

from .ba10c import viterbi
from .ba10e import print_mat
from .ba10h import estimate_tmat, estimate_emat


def parse_input(handle):
    niter = int(next(handle).rstrip())
    next(handle)
    seq = next(handle).rstrip()
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    states = next(handle).split()
    next(handle)
    lines = [next(handle) for _ in range(len(states) + 1)]
    tmat = {
        (states[i], states[j]): float(v)
        for i, x in enumerate(lines[1:])
        for j, v in enumerate(x.split()[1:])
    }
    next(handle)
    lines = [next(handle) for i in range(len(states) + 1)]
    emat = {
        (states[i], alphabet[j]): float(v)
        for i, x in enumerate(lines[1:])
        for j, v in enumerate(x.split()[1:])
    }
    return niter, seq, states, alphabet, tmat, emat


def main(file):
    niter, seq, st, al, tmat, emat = parse_input(open(file))
    for _ in range(niter):
        path = viterbi(seq, st, tmat, emat)
        tmat = estimate_tmat(path, st, to_dict=True)
        emat = estimate_emat(seq, al, path, st, to_dict=True)
    print_mat(estimate_tmat(path, st), st, st)
    print("--------")
    print_mat(estimate_emat(seq, al, path, st), st, al)
