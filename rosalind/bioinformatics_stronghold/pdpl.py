# Creating a Restriction Map

from .helpers import Parser

# The logic here is that, if we start at zero, the maximum of the differences
# observed in a difference multiset (ms), must exist as a number in the set x
# if we are starting at 0. Equally, all numbers of x must exist somewhere in the
# multiset (as differences from 0). So we can iterate through the multiset to
# find them.
#
# Since difference between all pairs of values in x must exist in the multiset,
# we can find a new value for x (from the multiset) by looking to see whether
# all differences between the new value (n) and all other known x values exist
# in the multiset.


def pdpl(ms):
    """Creating a Restriction Map"""
    x = [0, max(ms)]
    for n in set(ms):
        diffs = [abs(n - member) for member in x]
        if set(diffs).issubset(set(ms)):
            x.append(n)
            for d in diffs:
                ms.remove(d)

    return sorted(x)


def main(file):
    print(*pdpl(Parser(file).ints()))
