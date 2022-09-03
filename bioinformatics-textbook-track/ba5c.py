# Find a Longest Common Subsequence of Two Strings


def longest_common_subsequence(s1, s2):
    # Initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = 0
        p[j, 0] = "↑"

    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = "←"

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                opt = [m[j + 1, i], m[j, i + 1], m[j, i] + 1]
            else:
                opt = [m[j + 1, i], m[j, i + 1], m[j, i]]
            m[j + 1, i + 1] = max(opt)
            p[j + 1, i + 1] = ["←", "↑", "↖"][opt.index(max(opt))]

    # traceback
    ss = ""
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        if p[j, i] == "↖":
            ss += s1[i - 1]
            j, i = j - 1, i - 1
        elif p[j, i] == "←":
            i = i - 1
        elif p[j, i] == "↑":
            j = j - 1

    return ss[::-1]


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(longest_common_subsequence(s1, s2))
