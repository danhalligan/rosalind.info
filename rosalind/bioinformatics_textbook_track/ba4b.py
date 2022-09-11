# Compute the Number of Times a Pattern Appears in a Text

import re
from .ba4a import translate
from .ba1c import revcomp


def transcribe(dna):
    return re.sub("T", "U", dna)


def rev_transcribe(rna):
    return re.sub("U", "T", rna)


def find_matches(rna, pattern):
    for i in range(3):
        f = translate(rna[i:])
        for m in re.finditer(rf"(?=({pattern}))", f):
            start = m.span()[0] * 3 + i
            end = start + len(pattern) * 3
            yield rna[start:end]


def find_genome_substrings(dna, aa):
    for match in find_matches(transcribe(dna), aa):
        yield rev_transcribe(match)
    for match in find_matches(transcribe(revcomp(dna)), aa):
        yield revcomp(rev_transcribe(match))


def main(file):
    dna, aa = open(file).read().splitlines()
    for seq in find_genome_substrings(dna, aa):
        print(seq)
