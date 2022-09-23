# Implement PSMSearch

from .ba11f import peptide_identification


def main(file):
    *svs, proteome, T = open(file).read().splitlines()
    T = int(T)
    for sv in svs:
        sv = list(map(int, sv.split()))
        peptide, score = peptide_identification(sv, proteome)
        if score > T:
            print(peptide)
