# Implement DecodingIdealSpectrum

from .ba11a import spectrum_graph


def full(ions):
    def infer_peptide(w, seq, seen):
        if len(seq) == len(ions) // 2:
            yield seq
        for k in graph[w]:
            if k["n"] not in seen:
                yield from infer_peptide(k["n"], seq + k["l"], seen + [w])

    graph = spectrum_graph(ions)
    yield from infer_peptide(min(ions), "", [])


def main(file):
    x = [0] + list(map(int, open(file).read().split()))
    specs = list(full(x))
    print(specs[0])
