# Compute the Number of Peptides of Given Total Mass

from functools import cache
from .ba4c import mass

# The total number of ways to create a peptide of mass e.g. 1024 can be
# calculated as the sum over ways to create a peptide of mass (1024-aa) for each
# aa in the list of amino acid masses. These sub masses then be calculated
# recursively...

# Note that for the possible amino acid matches, we use the *set* of masses (not
# all masses, which differs because some amino acids have the same mass). This
# effectively means we do not count multiple solutions with different amino
# acids of the same weight! This is IMO a bug in rosalind as these should be
# counted.


def n_peptides(m):
    @cache
    def nways(m):
        if m < 0:
            return 0
        if m == 0:
            return 1
        return sum(nways(m - k) for k in masses)

    masses = set(mass().values())
    return nways(m)


def main(file):
    print(n_peptides(int(open(file).read())))
