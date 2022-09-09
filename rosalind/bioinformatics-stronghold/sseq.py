# Finding a Spliced Motif

from .helpers import Parser


def matches(s1, s2):
    i, j = 0, 0
    while j < len(s2):
        if s2[j] == s1[i]:
            yield i + 1
            j += 1
        i += 1


def main(file):
    s1, s2 = [x.seq for x in Parser(file).fastas()]
    print(*list(matches(s1, s2)))
