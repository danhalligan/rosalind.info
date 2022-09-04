# Implement ConvolutionCyclopeptideSequencing

from .ba4h import spectrum_convolution
from .ba4g import leaderboard_cyclopeptide_sequencing


def convolution_cyclopeptide_sequencing(m, n, spectrum):
    conv = spectrum_convolution(spectrum)
    conv = [(k, v) for k, v in conv if 57 <= k <= 200]
    masses = [k for k, v in conv if v >= conv[m - 1][1]]
    return leaderboard_cyclopeptide_sequencing(spectrum, n, masses)


def main(file):
    m, n, spectrum = open(file).read().splitlines()
    spectrum = list(map(int, spectrum.split()))
    p = convolution_cyclopeptide_sequencing(int(m), int(n), spectrum)
    print(*p, sep="-")
