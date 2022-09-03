# Find a Highest-Scoring Fitting Alignment of Two Strings


def fitting_alignment(s1, s2):
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = -j
        p[j, 0] = "↑"
    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = "←"

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            match = 1 if s1[i] == s2[j] else -1
            opt = [
                m[j, i] + match,
                m[j, i + 1] - 1,
                m[j + 1, i] - 1,
            ]
            m[new] = max(opt)
            p[new] = ["↖", "↑", "←"][opt.index(max(opt))]

    sc = [m[len(s2), i] for i in range(len(s1) + 1)]
    max_score = max(sc)
    i = sc.index(max_score)
    j = len(s2)

    a1, a2 = "", ""
    while i > 0 and j > 0:
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

    return max_score, a1[::-1], a2[::-1]


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(*fitting_alignment(s1, s2), sep="\n")
