# Genome Assembly as Shortest Superstring

from math import floor
from .helpers import Parser


# Fast find match using `find`
# Find overlap of at least min_overlap of prefix of s2 in s1
# Stores an index to search for next match
def find_overlap(s1, s2, min_overlap=None):
    ix = 1
    if min_overlap is None:
        min_overlap = floor(len(s2) / 2)
    while ix < len(s1):
        ix = s1.find(s2[:min_overlap], ix)
        if ix == -1:
            break
        if s2.startswith(s1[ix:]):
            return len(s1) - ix
        ix += 1


def construct_assembly(seqs):
    # Build a forward and reverse overlap graph of sequences
    fmap, rmap, starts, ends = ({}, {}, {}, {})
    for p1 in seqs.keys():
        for p2 in seqs.keys():
            if p1 in starts or p2 in ends or p1 in p2:
                continue
            n = find_overlap(seqs[p1], seqs[p2])
            if n:
                fmap[p1] = {"overlap": n, "next": p2}
                rmap[p2] = p1
                starts[p1] = True
                ends[p2] = True
                break

    # Find starting key using rmap
    k = list(seqs.keys())[0]
    while k in rmap:
        k = rmap[k]

    # Initialise list with sequence 1, then add suffixes of matching sequences
    seq = [seqs[k]]
    while k in fmap:
        seq.append(seqs[fmap[k]["next"]][fmap[k]["overlap"] :])
        k = fmap[k]["next"]

    # Join all sequences
    return "".join(seq)


def main(file):
    seqs = Parser(file).fastas()
    seqs = dict([(x.id, x.seq) for x in seqs])
    print(construct_assembly(seqs))
