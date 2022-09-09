# Edit Distance

from .helpers import Parser


def edit(s1, s2):
    """Edit Distance"""
    # initialise
    m = {}
    for j in range(len(s2) + 1):
        m[j, 0] = j
    for i in range(len(s1) + 1):
        m[0, i] = i

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
            else:
                m[j + 1, i + 1] = min([m[j + 1, i], m[j, i], m[j, i + 1]]) + 1

    return m[len(s2), len(s1)]


def main(file):
    seqs = Parser(file).seqs()
    print(edit(seqs[0], seqs[1]))
