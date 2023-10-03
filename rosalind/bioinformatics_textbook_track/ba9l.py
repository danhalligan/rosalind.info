# Implement BWMatching

from .ba9j import index_seq


def rindex(x, value):
    x.reverse()
    i = x.index(value)
    x.reverse()
    return len(x) - i - 1


def bwmatching(FirstColumn, LastColumn, Pattern):
    first = list(index_seq(FirstColumn))
    last = list(index_seq(LastColumn))

    top = 0
    bottom = len(LastColumn) - 1
    while top <= bottom:
        if Pattern:
            Pattern, symbol = Pattern[:-1], Pattern[-1]
            x = list(LastColumn[top : bottom + 1])
            if symbol in x:
                topIndex = x.index(symbol) + top
                bottomIndex = rindex(x, symbol) + top
                top = first.index(last[topIndex])
                bottom = first.index(last[bottomIndex])
            else:
                return 0
        else:
            return bottom - top + 1


def main(file):
    text, patterns = open(file).read().splitlines()
    patterns = patterns.split()
    print(*[bwmatching(sorted(text), text, pattern) for pattern in patterns])
