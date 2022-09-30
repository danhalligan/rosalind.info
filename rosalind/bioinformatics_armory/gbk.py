# GenBank Introduction

from Bio import Entrez
import os


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
    return Entrez.read(handle)


def main(file):
    orgm, start, end = open(file).read().splitlines()
    print(gbk(orgm, start, end)["Count"])
