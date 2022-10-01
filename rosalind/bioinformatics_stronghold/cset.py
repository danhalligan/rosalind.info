# Fixing an Inconsistent Character Set

# We can find conflicts if there's an intersection in all implied sets
# between two character rows

from collections import defaultdict
from itertools import product


def conflict(c1, c2):
    for a, b in product([1, 0], repeat=2):
        s1 = set(i for i, c in enumerate(c1) if c == a)
        s2 = set(i for i, c in enumerate(c2) if c == b)
        if len(s1.intersection(s2)) == 0:
            return False
    return True


def conflicts(characters):
    count = defaultdict(int)
    for i in range(len(characters)):
        for j in range(i + 1, len(characters)):
            if conflict(characters[i], characters[j]):
                count[i] += 1
                count[j] += 1
    return count


def main(file):
    lines = open(file).read().splitlines()
    characters = [[int(x) for x in list(ch)] for ch in lines]
    count = conflicts(characters)
    rm = [k for k, v in count.items() if v == max(count.values())][0]
    print(*(lines[:rm] + lines[rm + 1 :]), sep="\n")
