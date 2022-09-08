# Computing GC Content

from .helpers import Parser, Dna


def max_gc(seqs):
    gc = [Dna(rec.seq).gc_content() for rec in seqs]
    m = gc.index(max(gc))
    return {"name": seqs[m].id, "value": gc[m] * 100}


def main(file):
    res = max_gc(Parser(file).fastas())
    print(res["name"], round(res["value"], 5), sep="\n")
