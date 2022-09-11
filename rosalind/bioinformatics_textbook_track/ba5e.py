# Find a Highest-Scoring Alignment of Two Strings

from importlib import resources


def blosum62():
    path = resources.files("rosalind.resources").joinpath("blosum62.txt")
    lines = open(path).read().splitlines()
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])


def global_alignment(s1, s2, penalty=-5):
    score = blosum62()
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = penalty * j
        p[j, 0] = "↑"
    for i in range(len(s1) + 1):
        m[0, i] = penalty * i
        p[0, i] = "←"

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            opt = [
                m[j, i] + score[s1[i]][s2[j]],
                m[j, i + 1] + penalty,
                m[j + 1, i] + penalty,
            ]
            m[new] = max(opt)
            p[new] = ["↖", "↑", "←"][opt.index(max(opt))]

    i, j = len(s1), len(s2)
    a1, a2 = "", ""
    while i > 0 or j > 0:
        if p[j, i] == "↖":
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j, i = j - 1, i - 1
        elif p[j, i] == "←":
            a1 += s1[i - 1]
            a2 += "-"
            i = i - 1
        elif p[j, i] == "↑":
            a1 += "-"
            a2 += s2[j - 1]
            j = j - 1

    return m[len(s2), len(s1)], a1[::-1], a2[::-1]


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(*global_alignment(s1, s2), sep="\n")
