# Construct a Profile HMM with Pseudocounts

import numpy as np
from .ba10e import normalise, print_tprob, print_eprob, profile_hmm


def parse_input(handle):
    θ, sig = map(float, next(handle).rstrip().split())
    next(handle)
    alphabet = next(handle).split()
    next(handle)
    alignment = np.array([list(x) for x in handle.read().splitlines()])
    return θ, sig, alphabet, alignment


# Add pseudocounts to a transition probability matrix
# S  can go to I0, M1, D1
# In/Mn/Dn can go to In  M(n+1)  D(n+1)
# The indexing here is horrible...
def add_transition_pseudocounts(x, sig):
    n = (x.shape[0] - 3) // 3
    x[0, 1:4] += sig
    x[1, 1:4] += sig
    for i in range(0, n):
        x[i * 3 + 2 : i * 3 + 5, (i + 1) * 3 + 1 : (i + 1) * 3 + 4] += sig
    return normalise(x)


# Add pseudocounts to an emission probability matrix
def add_emission_pseudocounts(x, sig):
    n = (x.shape[0] - 3) // 3
    x[1, :] += sig
    for i in range(0, n):
        x[i * 3 + 2, :] += sig
        x[i * 3 + 4, :] += sig
    return normalise(x)


def main(file):
    θ, sig, alphabet, alignment = parse_input(open(file))
    tprob, eprob = profile_hmm(θ, alphabet, alignment)
    tprob = add_transition_pseudocounts(tprob, sig)
    eprob = add_emission_pseudocounts(eprob, sig)
    print_tprob(tprob)
    print("--------")
    print_eprob(eprob, alphabet)
