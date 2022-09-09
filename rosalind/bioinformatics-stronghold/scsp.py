# Interleaving Two Motifs

from .helpers import Parser


# We can calculate a shortest common supersequence by constructing a
# 2D matrix as before, but disallowing mismatches and with a special
# traceback.
def scsp(s1, s2):
    """Interleaving Two Motifs"""

    # Initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = j
        p[j, 0] = [j - 1, 0]

    for i in range(len(s1) + 1):
        m[0, i] = i
        p[0, i] = [0, i - 1]

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
                p[j + 1, i + 1] = [j, i]
            else:
                opt = [m[j + 1, i], m[j, i + 1]]
                m[j + 1, i + 1] = min(opt) + 1
                p[j + 1, i + 1] = [[j + 1, i], [j, i + 1]][opt.index(min(opt))]

    # traceback
    ss = ""
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        if p[j, i] == [j - 1, i - 1]:
            ss += s1[i - 1]
        elif p[j, i] == [j, i - 1]:
            ss += s1[i - 1]
        elif p[j, i] == [j - 1, i]:
            ss += s2[j - 1]
        j, i = p[j, i]

    return {"dist": m[len(s2), len(s1)], "ss": ss[::-1]}


def main(file):
    seqs = Parser(file).lines()
    out = scsp(seqs[0], seqs[1])
    print(out["ss"], sep="\n")
