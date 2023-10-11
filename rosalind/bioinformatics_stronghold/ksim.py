# Finding All Similar Motifs

# Here we use the fact that if two strings of length n match with at most d
# mismatches, then they must share a k-mer of length k = [n/(d + 1)]
# We can therefore take kmers from our pattern and search for perfect matches
# in the sequence.
#
# Since the pattern is long(ish) and the number of mismatches is relatively
# low, we can guarantee reasonably long perfect matching k-mers. These
# (potentially overlapping) perfect matches can be found quickly (with regex).
# In practice we have ~100 kmers of ~100bp for this problem
#
# Having found a seed, you might think we could extend them to find the maximum
# perfect matching section (and remove duplicates we find at his stage), but,
# we can't do this, as if we do, we can miss some non-optimal matches...
# Instead we need to extend our seeds incorporating mismatches / indels from
# the start and end.
#
# To approximate extensions in both directions, we:
# - consider either adding a match, an insertion in s or insertion in t
# - recurse for each possibility
# - stop recursion when reach end of sequence or score > k
# - combine start and end extensions in all combinations where sum scores <= k

# Performance is a problem here and we implement two speed ups:
# 1. parallel computation
# 2. for a given extension, once we've seen a combination of position in
#    s, position in t and score once, we don't need to consider this again.

import sys
import re
import multiprocessing as mp


def get_seeds(x, seq, k):
    seed_size = len(x) // (k + 1)
    for s1 in range(0, len(x) - seed_size + 1, seed_size):
        px = (s1, s1 + seed_size)
        seed = x[px[0] : px[1]]
        for m in re.finditer(rf"(?=({seed}))", seq):
            ps = (m.span()[0], m.span()[0] + seed_size)
            yield (px, ps)


def process_seed(args):
    def extend_fwd(i, j, score):
        if (i, j, score) not in seen:
            seen.update([(i, j, score)])
            if score <= k:
                if i == len(x) - 1:
                    yield i, j, score
                if i + 1 < len(x):
                    yield from extend_fwd(i + 1, j, score + 1)
                if j + 1 < len(seq):
                    yield from extend_fwd(i, j + 1, score + 1)
                if i + 1 < len(x) and j + 1 < len(seq):
                    yield from extend_fwd(
                        i + 1, j + 1, score + int(x[i + 1] != seq[j + 1])
                    )

    def extend_rev(i, j, score):
        if (i, j, score) not in seen:
            seen.update([(i, j, score)])
            if score <= k:
                if i == 0:
                    yield i, j, score
                if i - 1 >= 0:
                    yield from extend_rev(i - 1, j, score + 1)
                if j - 1 >= 0:
                    yield from extend_rev(i, j - 1, score + 1)
                if i - 1 >= 0 and j - 1 >= 0:
                    yield from extend_rev(
                        i - 1, j - 1, score + int(x[i - 1] != seq[j - 1])
                    )

    print(".", end="", file=sys.stderr)
    sys.stderr.flush()
    sys.setrecursionlimit(10000)
    seed, k, x, seq = args
    xcoord, seqcoord = seed
    res = set()
    seen = set()
    fwds = list(extend_fwd(xcoord[1] - 1, seqcoord[1] - 1, 0))
    if not fwds:
        return set()
    seen = set()
    revs = list(extend_rev(xcoord[0], seqcoord[0], 0))
    if not revs:
        return set()
    for i0, j0, s0 in revs:
        for i1, j1, s1 in fwds:
            if s0 + s1 <= k:
                res.add((j0 + 1, j1 - j0 + 1))
    return res


def main(file):
    k, x, seq = open(file).read().splitlines()
    k = int(k)
    seeds = list(get_seeds(x, seq, k))
    print(f"found {len(seeds)} seeds", file=sys.stderr)

    pool = mp.Pool(mp.cpu_count())
    args = ([seed, k, x, seq] for seed in seeds)
    res = pool.map(process_seed, args)
    res = set().union(*res)
    for x in sorted(list(res)):
        print(*x)
