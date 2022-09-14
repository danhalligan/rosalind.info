# Pattern Matching with the Suffix Array

from .ba9g import suffix_array


def pattern_matching_with_suffix_array(text, pattern, sa):
    minIndex = 0
    maxIndex = len(text)
    while minIndex < maxIndex:
        midIndex = (minIndex + maxIndex) // 2
        if pattern > text[sa[midIndex] :][: len(pattern)]:
            minIndex = midIndex + 1
        else:
            maxIndex = midIndex
    first = minIndex
    maxIndex = len(text)
    while minIndex < maxIndex:
        midIndex = (minIndex + maxIndex) // 2
        if pattern < text[sa[midIndex] :][: len(pattern)]:
            maxIndex = midIndex
        else:
            minIndex = midIndex + 1
    last = maxIndex
    if first > last:
        return []
    else:
        return list(range(first, last))


def main(file):
    seq, *patterns = open(file).read().splitlines()
    sa = suffix_array(seq)
    inds = []
    for p in patterns:
        for ind in pattern_matching_with_suffix_array(seq, p, sa):
            inds += [sa[ind]]
    print(*sorted(set(inds)))
