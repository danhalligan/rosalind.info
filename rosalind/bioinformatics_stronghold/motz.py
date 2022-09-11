# Motzkin Numbers and RNA Secondary Structures

from functools import cache
from .helpers import Parser
from .cat import valid_pair


@cache
def motz(seq, mod=10**6):
    """Motzkin Numbers and RNA Secondary Structures"""
    if len(seq) in range(1):
        return 1
    else:
        return (
            motz(seq[1:])
            + sum(
                motz(seq[1:m]) * motz(seq[m + 1 :])
                for m in range(1, len(seq))
                if valid_pair(seq[0], seq[m])
            )
            % mod
        )


def main(file):
    print(motz(Parser(file).seqs()[0]))
