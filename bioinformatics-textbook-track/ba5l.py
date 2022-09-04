# Align Two Strings Using Linear Space

from .ba5e import blosum62
from .ba5k import middle_edge


def alignment_score(s1, s2, scores, penalty):
    return sum(
        penalty if s1[i] == "-" or s2[i] == "-" else scores[s1[i]][s2[i]]
        for i in range(len(s1))
    )


def find_path(s1, s2, scores, penalty):
    def lsa(t, b, left, right):
        if left == right:
            return "↓" * (b - t)
        elif t == b:
            return "→" * (right - left)
        else:
            ((i, j), (i2, j2)) = middle_edge(s1[t:b], s2[left:right], scores, penalty)
            edge = "↓" if j == j2 else "→" if i == i2 else "↘"
            return (
                lsa(t, i + t, left, j + left) + edge + lsa(i2 + t, b, j2 + left, right)
            )

    return lsa(0, len(s1), 0, len(s2))


def build_alignment(path, s1, s2):
    a1, a2 = "", ""
    i, j = 0, 0
    for arrow in path:
        if arrow == "↘":
            a1 += s1[i]
            a2 += s2[j]
            i += 1
            j += 1
        elif arrow == "↓":
            a1 += s1[i]
            a2 += "-"
            i += 1
        else:
            a1 += "-"
            a2 += s2[j]
            j += 1
    return a1, a2


def main(file):
    s1, s2 = open(file).read().splitlines()
    scores = blosum62()
    path = find_path(s1, s2, scores, -5)
    a1, a2 = build_alignment(path, s1, s2)
    print(alignment_score(a1, a2, scores, -5))
    print(a1, a2, sep="\n")
