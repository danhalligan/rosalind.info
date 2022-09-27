# Implement Baum-Welch Learning

import numpy as np
from .ba10i import parse_input, print_dict
from .ba10j import soft_decode, forward, backward
from .ba10h import as_dict

# To implement this method, we need to calculate two responsibility matrices
# (PI^* and PI^**, see page 225). PI^* is the same as `soft_decode` calculated
# in ba10j.


# Estimate responsibility matrix (PI^**)
# implementing equation in purple box on p. 224
# NB: Weight_i(l,k) = transition_pi_(i-1),pi_i * emission_pi_i(x_i)
def estimate_pi2(seq, fwd, bak, tmat, emat, states):
    rep_mat = np.zeros((fwd.shape[0] - 1, len(states), len(states)), dtype=float)
    for i in range(0, fwd.shape[0] - 1):
        for j, s1 in enumerate(states):
            for k, s2 in enumerate(states):
                weight = tmat[s1, s2] * emat[s2, seq[i + 1]]
                rep_mat[i, j, k] = (
                    fwd[i, j] * bak[i + 1, k] * weight / sum(fwd[i, :] * bak[i, :])
                )
    return rep_mat


# Estimate a transition matrix from the responsibility matrix PI^**
def estimate_tmat(seq, st, tmat, emat):
    fwd = forward(seq, st, tmat, emat)
    bak = backward(seq, st, tmat, emat)
    pi2 = estimate_pi2(seq, fwd, bak, tmat, emat, st)
    tmat = np.sum(pi2, 0)
    tmat = tmat / np.sum(tmat, axis=1, keepdims=True)
    return as_dict(tmat, st, st)


# Estimate a emission matrix from the responsibility matrix PI^*
def estimate_emat(seq, al, st, tmat, emat):
    pi1 = soft_decode(seq, st, tmat, emat)
    emat = np.zeros((len(st), len(al)), dtype=float)
    for i, emission in enumerate(al):
        ind = np.array(list(seq)) == emission
        emat[:, i] = np.sum(pi1[ind, :], 0)
    emat = emat / np.sum(emat, axis=1, keepdims=True)
    return as_dict(emat, st, al)


def main(file):
    niter, seq, st, al, tmat, emat = parse_input(open(file))
    for _ in range(niter):
        tmat2 = estimate_tmat(seq, st, tmat, emat)
        emat2 = estimate_emat(seq, al, st, tmat, emat)
        emat, tmat = emat2, tmat2
    print_dict(tmat, st, st)
    print("--------")
    print_dict(emat, st, al)
