# Reversal Distance

from itertools import combinations


def flip(x, i, j):
    """Flip a section of a sequence"""
    rev = list.copy(x)
    rev[i:j] = rev[i:j][::-1]
    return rev


def breaks(s, t):
    """Identify breaks between a sequence and target"""
    return [
        i + 1 for i in range(len(s) - 1) if abs(t.index(s[i]) - t.index(s[i + 1])) != 1
    ]


def min_breaks(seqs, t):
    rev = []
    for s in seqs:
        for i, j in combinations(breaks(s["s"], t), 2):
            rev.append({"s": flip(s["s"], i, j), "p": s["p"] + [[i, j - 1]]})
    min_b = len(t)
    mr = []
    for r in rev:
        n = len(breaks(r["s"], t))
        if n < min_b:
            min_b = n
            mr = [r]
        elif n == min_b:
            mr.append(r)
    return mr


# based on https://medium.com/@matthewwestmk/87c62d690eef
def sort(s, t):
    """Sorting by Reversals"""
    s = ["-"] + s + ["+"]
    t = ["-"] + t + ["+"]
    nr = 0
    c = [{"s": s, "p": []}]
    seqs = [s]
    while t not in seqs:
        c = min_breaks(c, t)
        nr += 1
        seqs = [x["s"] for x in c]
    return nr, c


def main(file):
    data = open(file).read().strip().split("\n\n")
    data = [tuple([list(map(int, y.split())) for y in x.split("\n")]) for x in data]
    print(*[sort(s, t)[0] for s, t in data])
