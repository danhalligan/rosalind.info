# Consensus and Profile

import numpy as np
from .helpers import Parser


def profile_matrix(seqs):
    def count_bases(v):
        return [sum(v == b) for b in "ACGT"]

    x = np.array([list(x) for x in seqs])
    return np.array([count_bases(v) for v in x.T]).T


def consensus_sequence(mat):
    return "".join(["ACGT"[np.argmax(v)] for v in mat.T])


def main(file):
    x = Parser(file).seqs()
    mat = profile_matrix(x)
    print(consensus_sequence(mat))
    for i in range(0, 4):
        print("ACGT"[i] + ":", *mat[i])
