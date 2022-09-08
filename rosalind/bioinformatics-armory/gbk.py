# GenBank Introduction

from Bio import Entrez
import os
import sys


def entrez_email():
    try:
        email = os.environ["ENTREZ_EMAIL"]
    except KeyError:
        print("Please add your email as the environment variable 'ENTREZ_EMAIL'")
        exit()
    return email


def gbk(orgm, start, end):
    Entrez.email = entrez_email()
    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{orgm}[Organism]",
        datetype="pdat",
        mindate=start,
        maxdate=end,
    )
    record = Entrez.read(handle)
    return record["Count"]


def main(file):
    orgm, start, end = open(file).read().splitlines()
    print(gbk(orgm, start, end))


if __name__ == "__main__":
    main(sys.argv[1])
