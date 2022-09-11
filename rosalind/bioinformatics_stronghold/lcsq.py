# Finding a Shared Spliced Motif

from .helpers import Parser


def lcsq(s1, s2):
    """Finding a Shared Spliced Motif"""
    # initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = 0
        p[j, 0] = [j - 1, 0]

    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = [0, i - 1]

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i] + 1
                p[j + 1, i + 1] = [j, i]
            else:
                opt = [m[j + 1, i], m[j, i + 1]]
                m[j + 1, i + 1] = max(opt)
                p[j + 1, i + 1] = [[j + 1, i], [j, i + 1]][opt.index(max(opt))]

    # traceback
    subs = ""
    i, j = len(s1), len(s2)
    while i > 0 and j > 0:
        if p[j, i] == [j - 1, i - 1]:
            subs += s1[i - 1]
        j, i = p[j, i]

    return subs[::-1]


def main(file):
    print(lcsq(*Parser(file).seqs()))
