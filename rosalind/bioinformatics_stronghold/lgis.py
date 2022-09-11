# Longest Increasing Subsequence

from .helpers import Parser


def lgis(x):
    """DP approach to longest increasing subsequence"""
    n = len(x)
    d = [1] * n  # length of the longest increasing subsequence ending at i
    p = [-1] * n  # pointer to previous part of subsequence

    for i in range(n):
        for j in range(i):
            if x[j] < x[i] and d[i] < d[j] + 1:
                d[i] = d[j] + 1
                p[i] = j

    ans = max(d)
    pos = d.index(ans)

    # traceback
    subseq = []
    while pos != -1:
        subseq.append(x[pos])
        pos = p[pos]

    subseq.reverse()
    return subseq


def main(file):
    data = Parser(file).lines()[1]
    data = [int(x) for x in data.split(" ")]
    s1 = lgis(data)
    print(*s1)

    s2 = lgis([-x for x in data])
    s2 = [-x for x in s2]
    print(*s2)
