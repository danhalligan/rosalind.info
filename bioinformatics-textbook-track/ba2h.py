# Implement DistanceBetweenPatternAndStrings

from .ba2b import distance_between_pattern_and_strings


def main(file):
    pattern, dna = open(file).read().splitlines()
    print(distance_between_pattern_and_strings(pattern, dna.split()))
