# Translate an RNA String into an Amino Acid String

import yaml
import re
import math
from importlib import resources


def genetic_code():
    path = resources.files("rosalind.resources").joinpath("genetic_code.yaml")
    with open(path) as stream:
        return yaml.safe_load(stream)


def translate(rna):
    code = genetic_code()
    end = math.floor(len(rna) / 3) * 3
    prot = [code[rna[i : i + 3]] for i in range(0, end, 3)]
    return re.sub("\\*$", "", "".join(prot))


def main(file):
    rna = open(file).read().splitlines()[0]
    print(translate(rna))
