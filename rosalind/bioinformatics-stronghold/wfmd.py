# The Wright-Fisher Model of Genetic Drift

from .helpers import Parser, dbinom
import numpy as np


def wf_model(n, m, g):
    # transition matrix
    tmat = [[dbinom(x, n, m / n) for x in range(n + 1)] for m in range(n + 1)]
    tmat = np.array(tmat)
    v = np.array([0] * (n + 1))
    v[m] = 1
    for i in range(g):
        v = np.dot(v, tmat)
    return v


def wfmd(n, m, g, k):
    p = wf_model(2 * n, m, g)
    return sum(p[:-k])


def main(file):
    print(round(wfmd(*Parser(file).ints()), 3))
