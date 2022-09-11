# Finding a Protein Motif

import re
import requests as r
from io import StringIO
from .helpers import Parser, read_fasta


def find_protein_motif(seq, pattern="N[^P][ST][^P]"):
    motif = re.compile("(?=(" + pattern + "))")
    return [m.start() + 1 for m in motif.finditer(seq)]


def uniprot_output(id):
    return r.post(f"https://www.uniprot.org/uniprot/{id}.fasta").text


def get_uniprot(id):
    txt = uniprot_output(id)
    return list(read_fasta(StringIO(txt)))[0]


def main(file):
    for id in Parser(file).lines():
        seq = get_uniprot(id)
        matches = find_protein_motif(str(seq.seq))
        if len(matches):
            print(id)
            print(*matches)
