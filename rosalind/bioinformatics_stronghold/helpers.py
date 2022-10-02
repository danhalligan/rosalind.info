import re
import string
from rosalind.helpers import genetic_code
from math import comb


class Parser:
    """Parse problem data text files"""

    def __init__(self, file):
        if type(file) == str:
            self.handle = open(file)
        else:
            self.handle = file

    def dna(self):
        """Return the first line as a DNA string"""
        return Dna(self.line())

    def rna(self):
        """Return the first line as a RNA string"""
        return Rna(self.line())

    def line(self):
        """Return the first line"""
        return self.handle.readline().rstrip()

    def lines(self):
        """Return lines as a chomped list"""
        return self.handle.read().splitlines()

    def fastas(self):
        """Return fasta records as a list"""
        return list(read_fasta(self.handle))

    def seqs(self):
        """Return sequences from fasta records as a list"""
        return [x.seq for x in self.fastas()]

    def ints(self):
        """Return space separated integers from first line"""
        return list(map(int, self.line().split()))

    def floats(self):
        """Return space separated floats from first line"""
        return list(map(float, self.line().split()))


class Rec:
    """A simple FASTA record"""

    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class Seq:
    """A sequence with a specified alphabet (either DNA, RNA or Protein)"""

    def __init__(self, seq: str):
        if 0 in [c in self.alphabet for c in seq]:
            raise TypeError("String contains invalid characters!")
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __eq__(self, other):
        return self.seq == other

    def __getitem__(self, value):
        return type(self)(self.seq.__getitem__(value))

    def __repr__(self):
        return self.__class__.__name__ + "(" + self.seq + ")"

    def table(self):
        """Count number of each letter"""
        return {x: self.seq.count(x) for x in self.alphabet}

    @property
    def alphabet(self):
        return string.ascii_uppercase


class Dna(Seq):
    """A DNA sequence"""

    def rna(self):
        """Convert DNA to RNA"""
        return Rna(re.sub("T", "U", self.seq))

    def revc(self):
        """Reverse complement"""
        return Dna(self.seq[::-1].translate(str.maketrans("ACGT", "TGCA")))

    def gc_content(self):
        """Calculate GC content"""
        return sum([self.seq.count(base) for base in "GC"]) / len(self)

    def translate(self):
        """Translate DNA to protein sequence. The stop codon is removed
        automatically."""
        return self.rna().translate()

    @property
    def alphabet(self):
        return "ACGT"


class Rna(Seq):
    """An RNA sequence"""

    def __init__(self, seq: str):
        super().__init__(seq)
        self._code = self._genetic_code()

    @property
    def alphabet(self):
        return "ACGU"

    def translate(self):
        """Translate RNA to protein sequence. The stop codon is removed
        automatically."""
        x = self.seq
        prot = [self._code[x[i : i + 3]] for i in range(0, len(x), 3)]
        prot = "".join(prot)
        return Prot(re.sub("\\*$", "", prot))

    def _genetic_code(self):
        return genetic_code()


class Prot(Seq):
    """A protein sequence"""

    @property
    def alphabet(self):
        return "*ACDEFGHIKLMNPQRSTVWY"


def read_fasta(handle):
    header, sequence = "", []
    for line in handle:
        if line[0] == ">":
            if sequence:
                yield Rec(header, "".join(sequence))
            header, sequence = line[1:-1], []
        else:
            sequence.append(line.strip())
    yield Rec(header, "".join(sequence))


def dbinom(x, size, prob):
    """Binomial density"""
    return comb(size, x) * prob**x * (1.0 - prob) ** (size - x)


def pbinom(q, size, prob):
    """Binomial distribution function"""
    return sum([dbinom(x, size, prob) for x in range(0, q + 1)])


class Graph:
    """Simple graph implementation"""

    def __init__(self):
        self.g = {}

    def add_edge(self, a, b, weight=1, label=""):
        if a not in self.g:
            self.g[a] = {}
        self.g.append({"n": b, "w": weight, "l": label})

    def nodes(self):
        s = list(self.g.keys())
        e = [y["n"] for v in self.g.values() for y in v]
        return set(s) | set(e)


def nodes(graph):
    s = list(graph.keys())
    e = [y["n"] for v in graph.values() for y in v]
    return set(s) | set(e)
