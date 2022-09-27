![http://rosalind.info](/images/logo.png)

This repository contains solutions to bioinformatics coding challenges from
[rosalind.info]. Problems are organised by the various different
[locations]:

* [Python Village]: initial problems to learn a few basics about the Python
  programming language.
* [Bioinformatics Stronghold]: problems to discover the algorithms underlying a
  variety of bioinformatics topics.
* [Bioinformatics Armory]: unlike the stronghold in the Armory we solve
* [Bioinformatics Textbook Track]: problems associated with [Bioinformatics
  Algorithms: An Active Learning Approach].
* [Algorithmic Heights]: exercises to accompany the book [Algorithms].
  problems using existing tools.

## Running the solutions

This repository is written as a python module and uses [poetry] and [typer].

Solutions for each problem are located in individual files inside the directory
for each location.

You can install the versions of dependencies used here with:

```{shell}
poetry install
```

To run solutions within this environment run, e.g.:

```{shell}
poetry run rosalind ini2 rosalind_ini2.txt
```

To run the solution on the provided "Sample Dataset" from [rosalind.info] (which
should reproduce the "Sample Output"), run the solution in "test" mode:

```{shell}
poetry run rosalind --test ini2
```

## Testing

`pytest-snapshot` is used to test solutions to problems. In many cases solutions
generated will and should exactly match the "Sample Output" given at
rosalind.info. In cases cases, where e.g. ordering is not important, the
expected solutions (in `tests/expected`) have been updated to match code used
here, but are equally valid solutions.

To run the tests use:

```{shell}
poetry run pytest
```

To update the tests (adding or modifying snapshots / expected output) use:

```{shell}
poetry run pytest --snapshot-update
```

Note that some solutions (that use Entrez) require an email address. This
should be set as an environment variable, e.g.:

```{shell}
export ENTREZ_EMAIL=rosalind.franklin@cam.ac.uk
```

## About

* My rosalind profile: https://rosalind.info/users/danhalligan/

## Solutions

### Python Village

- [x] [INI1: Installing Python](rosalind/python_village/ini1.py)
- [x] [INI2: Variables and Some Arithmetic](rosalind/python_village/ini2.py)
- [x] [INI3: Strings and Lists](rosalind/python_village/ini3.py)
- [x] [INI4: Conditions and Loops](rosalind/python_village/ini4.py)
- [x] [INI5: Working with Files](rosalind/python_village/ini5.py)
- [x] [INI6: Dictionaries](rosalind/python_village/ini6.py)

### Bioinformatics Stronghold

- [x] [DNA: Counting DNA Nucleotides](rosalind/bioinformatics_stronghold/dna.py)
- [x] [RNA: Transcribing DNA into RNA](rosalind/bioinformatics_stronghold/rna.py)
- [x] [REVC: Complementing a Strand of DNA](rosalind/bioinformatics_stronghold/revc.py)
- [x] [FIB: Rabbits and Recurrence Relations](rosalind/bioinformatics_stronghold/fib.py)
- [x] [GC: Computing GC Content](rosalind/bioinformatics_stronghold/gc.py)
- [x] [HAMM: Counting Point Mutations](rosalind/bioinformatics_stronghold/hamm.py)
- [x] [IPRB: Mendel's First Law](rosalind/bioinformatics_stronghold/iprb.py)
- [x] [PROT: Translating RNA into Protein](rosalind/bioinformatics_stronghold/prot.py)
- [x] [SUBS: Finding a Motif in DNA](rosalind/bioinformatics_stronghold/subs.py)
- [x] [CONS: Consensus and Profile](rosalind/bioinformatics_stronghold/cons.py)
- [x] [FIBD: Mortal Fibonacci Rabbits](rosalind/bioinformatics_stronghold/fibd.py)
- [x] [GRPH: Overlap Graphs](rosalind/bioinformatics_stronghold/grph.py)
- [x] [IEV: Calculating Expected Offspring](rosalind/bioinformatics_stronghold/iev.py)
- [x] [LCSM: Finding a Shared Motif](rosalind/bioinformatics_stronghold/lcsm.py)
- [x] [LIA: Independent Alleles](rosalind/bioinformatics_stronghold/lia.py)
- [x] [MPRT: Finding a Protein Motif](rosalind/bioinformatics_stronghold/mprt.py)
- [x] [MRNA: Inferring mRNA from Protein](rosalind/bioinformatics_stronghold/mrna.py)
- [x] [ORF: Open Reading Frames](rosalind/bioinformatics_stronghold/orf.py)
- [x] [PERM: Enumerating Gene Orders](rosalind/bioinformatics_stronghold/perm.py)
- [x] [PRTM: Calculating Protein Mass](rosalind/bioinformatics_stronghold/prtm.py)
- [x] [REVP: Locating Restriction Sites](rosalind/bioinformatics_stronghold/revp.py)
- [x] [SPLC: RNA Splicing](rosalind/bioinformatics_stronghold/splc.py)
- [x] [LEXF: Enumerating k-mers Lexicographically](rosalind/bioinformatics_stronghold/lexf.py)
- [x] [LGIS: Longest Increasing Subsequence](rosalind/bioinformatics_stronghold/lgis.py)
- [x] [LONG: Genome Assembly as Shortest Superstring](rosalind/bioinformatics_stronghold/long.py)
- [x] [PMCH: Perfect Matchings and RNA Secondary Structures](rosalind/bioinformatics_stronghold/pmch.py)
- [x] [PPER: Partial Permutations](rosalind/bioinformatics_stronghold/pper.py)
- [x] [PROB: Introduction to Random Strings](rosalind/bioinformatics_stronghold/prob.py)
- [x] [SIGN: Enumerating Oriented Gene Orderings](rosalind/bioinformatics_stronghold/sign.py)
- [x] [SSEQ: Finding a Spliced Motif](rosalind/bioinformatics_stronghold/sseq.py)
- [x] [TRAN: Transitions and Transversions](rosalind/bioinformatics_stronghold/tran.py)
- [x] [TREE: Completing a Tree](rosalind/bioinformatics_stronghold/tree.py)
- [x] [CAT: Catalan Numbers and RNA Secondary Structures](rosalind/bioinformatics_stronghold/cat.py)
- [x] [CORR: Error Correction in Reads](rosalind/bioinformatics_stronghold/corr.py)
- [x] [INOD: Counting Phylogenetic Ancestors](rosalind/bioinformatics_stronghold/inod.py)
- [x] [KMER: k-Mer Composition](rosalind/bioinformatics_stronghold/kmer.py)
- [x] [KMP: Speeding Up Motif Finding](rosalind/bioinformatics_stronghold/kmp.py)
- [x] [LCSQ: Finding a Shared Spliced Motif](rosalind/bioinformatics_stronghold/lcsq.py)
- [x] [LEXV: Ordering Strings of Varying Length Lexicographically](rosalind/bioinformatics_stronghold/lexv.py)
- [x] [MMCH: Maximum Matchings and RNA Secondary Structures](rosalind/bioinformatics_stronghold/mmch.py)
- [x] [PDST: Creating a Distance Matrix](rosalind/bioinformatics_stronghold/pdst.py)
- [x] [REAR: Reversal Distance](rosalind/bioinformatics_stronghold/rear.py)
- [x] [RSTR: Matching Random Motifs](rosalind/bioinformatics_stronghold/rstr.py)
- [x] [SSET: Counting Subsets](rosalind/bioinformatics_stronghold/sset.py)
- [x] [ASPC: Introduction to Alternative Splicing](rosalind/bioinformatics_stronghold/aspc.py)
- [x] [EDIT: Edit Distance](rosalind/bioinformatics_stronghold/edit.py)
- [x] [EVAL: Expected Number of Restriction Sites](rosalind/bioinformatics_stronghold/eval.py)
- [x] [MOTZ: Motzkin Numbers and RNA Secondary Structures](rosalind/bioinformatics_stronghold/motz.py)
- [x] [SCSP: Interleaving Two Motifs](rosalind/bioinformatics_stronghold/scsp.py)
- [x] [SETO: Introduction to Set Operations](rosalind/bioinformatics_stronghold/seto.py)
- [x] [SORT: Sorting by Reversals](rosalind/bioinformatics_stronghold/sort.py)
- [x] [SPEC: Inferring Protein from Spectrum](rosalind/bioinformatics_stronghold/spec.py)
- [x] [TRIE: Introduction to Pattern Matching](rosalind/bioinformatics_stronghold/trie.py)
- [x] [CONV: Comparing Spectra with the Spectral Convolution](rosalind/bioinformatics_stronghold/conv.py)
- [x] [DBRU: Constructing a De Bruijn Graph](rosalind/bioinformatics_stronghold/dbru.py)
- [x] [EDTA: Edit Distance Alignment](rosalind/bioinformatics_stronghold/edta.py)
- [x] [FULL: Inferring Peptide from Full Spectrum](rosalind/bioinformatics_stronghold/full.py)
- [x] [INDC: Independent Segregation of Chromosomes](rosalind/bioinformatics_stronghold/indc.py)
- [x] [LREP: Finding the Longest Multiple Repeat](rosalind/bioinformatics_stronghold/lrep.py)
- [x] [RNAS: Wobble Bonding and RNA Secondary Structures](rosalind/bioinformatics_stronghold/rnas.py)
- [x] [AFRQ: Counting Disease Carriers](rosalind/bioinformatics_stronghold/afrq.py)
- [x] [CTEA: Counting Optimal Alignments](rosalind/bioinformatics_stronghold/ctea.py)
- [x] [GLOB: Global Alignment with Scoring Matrix](rosalind/bioinformatics_stronghold/glob.py)
- [x] [PCOV: Genome Assembly with Perfect Coverage](rosalind/bioinformatics_stronghold/pcov.py)
- [x] [PRSM: Matching a Spectrum to a Protein](rosalind/bioinformatics_stronghold/prsm.py)
- [x] [SGRA: Using the Spectrum Graph to Infer Peptides](rosalind/bioinformatics_stronghold/sgra.py)
- [x] [SUFF: Encoding Suffix Trees](rosalind/bioinformatics_stronghold/suff.py)
- [x] [GASM: Genome Assembly Using Reads](rosalind/bioinformatics_stronghold/gasm.py)
- [x] [GCON: Global Alignment with Constant Gap Penalty](rosalind/bioinformatics_stronghold/gcon.py)
- [x] [LING: Linguistic Complexity of a Genome](rosalind/bioinformatics_stronghold/ling.py)
- [x] [LOCA: Local Alignment with Scoring Matrix](rosalind/bioinformatics_stronghold/loca.py)
- [x] [MGAP: Maximizing the Gap Symbols of an Optimal Alignment](rosalind/bioinformatics_stronghold/mgap.py)
- [x] [MULT: Multiple Alignment](rosalind/bioinformatics_stronghold/mult.py)
- [x] [PDPL: Creating a Restriction Map](rosalind/bioinformatics_stronghold/pdpl.py)
- [x] [SEXL: Sex-Linked Inheritance](rosalind/bioinformatics_stronghold/sexl.py)
- [x] [WFMD: The Wright-Fisher Model of Genetic Drift](rosalind/bioinformatics_stronghold/wfmd.py)
- [x] [ASMQ: Assessing Assembly Quality with N50 and N75](rosalind/bioinformatics_stronghold/asmq.py)
- [x] [EBIN: Wright-Fisher's Expected Behavior](rosalind/bioinformatics_stronghold/ebin.py)
- [x] [FOUN: The Founder Effect and Genetic Drift](rosalind/bioinformatics_stronghold/foun.py)
- [x] [GAFF: Global Alignment with Scoring Matrix and Affine Gap Penalty](rosalind/bioinformatics_stronghold/gaff.py)
- [x] [GREP: Genome Assembly with Perfect Coverage and Repeats](rosalind/bioinformatics_stronghold/grep.py)
- [x] [OAP: Overlap Alignment](rosalind/bioinformatics_stronghold/oap.py)
- [x] [SIMS: Finding a Motif with Modifications](rosalind/bioinformatics_stronghold/sims.py)
- [x] [SMGB: Semiglobal Alignment](rosalind/bioinformatics_stronghold/smgb.py)
- [ ] NWCK: Distances in Trees
- [ ] ITWV: Finding Disjoint Motifs in a Gene
- [ ] MREP: Identifying Maximal Repeats
- [ ] KSIM: Finding All Similar Motifs
- [ ] LAFF: Local Alignment with Affine Gap Penalty
- [ ] OSYM: Isolating Symbols in Alignments
- [ ] CTBL: Creating a Character Table
- [ ] NKEW: Newick Format with Edge Weights
- [ ] CSTR: Creating a Character Table from Genetic Strings
- [ ] CUNR: Counting Unrooted Binary Trees
- [ ] QRT: Quartets
- [ ] CHBP: Character-Based Phylogeny
- [ ] CNTQ: Counting Quartets
- [ ] EUBT: Enumerating Unrooted Binary Trees
- [ ] MEND: Inferring Genotype from a Pedigree
- [ ] ROOT: Counting Rooted Binary Trees
- [ ] SPTD: Phylogeny Comparison with Split Distance
- [ ] ALPH: Alignment-Based Phylogeny
- [ ] CSET: Fixing an Inconsistent Character Set
- [ ] QRTD: Quartet Distance
- [ ] RSUB: Identifying Reversing Substitutions

### Bioinformatics Armory

- [x] [INI: Introduction to the Bioinformatics Armory](rosalind/bioinformatics_armory/ini.py)
- [x] [GBK: GenBank Introduction](rosalind/bioinformatics_armory/gbk.py)
- [x] [FRMT: Data Formats](rosalind/bioinformatics_armory/frmt.py)
- [x] MEME: New Motif Discovery
- [x] [NEED: Pairwise Global Alignment](rosalind/bioinformatics_armory/need.py)
- [x] [TFSQ: FASTQ format introduction](rosalind/bioinformatics_armory/tfsq.py)
- [x] [PHRE: Read Quality Distribution](rosalind/bioinformatics_armory/phre.py)
- [x] [PTRA: Protein Translation](rosalind/bioinformatics_armory/ptra.py)
- [x] [FILT: Read Filtration by Quality](rosalind/bioinformatics_armory/filt.py)
- [x] [RVCO: Complementing a Strand of DNA](rosalind/bioinformatics_armory/rvco.py)
- [x] [SUBO: Suboptimal Local Alignment](rosalind/bioinformatics_armory/subo.py)
- [x] [BPHR: Base Quality Distribution](rosalind/bioinformatics_armory/bphr.py)
- [x] CLUS: Global Multiple Alignment
- [x] [ORFR: Finding Genes with ORFs](rosalind/bioinformatics_armory/orfr.py)
- [x] [BFIL: Base Filtration by Quality](rosalind/bioinformatics_armory/bfil.py)

#### Notes

* For "MEME" and "CLUS" I have not written a solution. Use the web interface as
  instructed.
* For "SUBO", you need to run the online interface, identify the 32-40 bp and
  then can use the solution here to count the occurrences of this in
  the sequences.

### Bioinformatics Textbook Track

- [x] [BA1A: Compute the Number of Times a Pattern Appears in a Text](rosalind/bioinformatics_textbook_track/ba1a.py)
- [x] [BA1B: Find the Most Frequent Words in a String](rosalind/bioinformatics_textbook_track/ba1b.py)
- [x] [BA1C: Find the Reverse Complement of a String](rosalind/bioinformatics_textbook_track/ba1c.py)
- [x] [BA1D: Find All Occurrences of a Pattern in a String](rosalind/bioinformatics_textbook_track/ba1d.py)
- [x] [BA1E: Find Patterns Forming Clumps in a String](rosalind/bioinformatics_textbook_track/ba1e.py)
- [x] [BA1F: Find a Position in a Genome Minimizing the Skew](rosalind/bioinformatics_textbook_track/ba1f.py)
- [x] [BA1G: Compute the Hamming Distance Between Two Strings](rosalind/bioinformatics_textbook_track/ba1g.py)
- [x] [BA1H: Find All Approximate Occurrences of a Pattern in a String](rosalind/bioinformatics_textbook_track/ba1h.py)
- [x] [BA1I: Find the Most Frequent Words with Mismatches in a String](rosalind/bioinformatics_textbook_track/ba1i.py)
- [x] [BA1J: Find Frequent Words with Mismatches and Reverse Complements](rosalind/bioinformatics_textbook_track/ba1j.py)
- [x] [BA1K: Generate the Frequency Array of a String](rosalind/bioinformatics_textbook_track/ba1k.py)
- [x] [BA1L: Implement PatternToNumber](rosalind/bioinformatics_textbook_track/ba1l.py)
- [x] [BA1M: Implement NumberToPattern](rosalind/bioinformatics_textbook_track/ba1m.py)
- [x] [BA1N: Generate the d-Neighborhood of a String](rosalind/bioinformatics_textbook_track/ba1n.py)
- [x] [BA2A: Implement MotifEnumeration](rosalind/bioinformatics_textbook_track/ba2a.py)
- [x] [BA2B: Find a Median String](rosalind/bioinformatics_textbook_track/ba2b.py)
- [x] [BA2C: Find a Profile-most Probable k-mer in a String](rosalind/bioinformatics_textbook_track/ba2c.py)
- [x] [BA2D: Implement GreedyMotifSearch](rosalind/bioinformatics_textbook_track/ba2d.py)
- [x] [BA2E: Implement GreedyMotifSearch with Pseudocounts](rosalind/bioinformatics_textbook_track/ba2e.py)
- [x] [BA2F: Implement RandomizedMotifSearch](rosalind/bioinformatics_textbook_track/ba2f.py)
- [x] [BA2G: Implement GibbsSampler](rosalind/bioinformatics_textbook_track/ba2g.py)
- [x] [BA2H: Implement DistanceBetweenPatternAndStrings](rosalind/bioinformatics_textbook_track/ba2h.py)
- [x] [BA3A: Generate the k-mer Composition of a String](rosalind/bioinformatics_textbook_track/ba3a.py)
- [x] [BA3B: Reconstruct a String from its Genome Path](rosalind/bioinformatics_textbook_track/ba3b.py)
- [x] [BA3C: Construct the Overlap Graph of a Collection of k-mers](rosalind/bioinformatics_textbook_track/ba3c.py)
- [x] [BA3D: Construct the De Bruijn Graph of a String](rosalind/bioinformatics_textbook_track/ba3d.py)
- [x] [BA3E: Construct the De Bruijn Graph of a Collection of k-mers](rosalind/bioinformatics_textbook_track/ba3e.py)
- [x] [BA3F: Find an Eulerian Cycle in a Graph](rosalind/bioinformatics_textbook_track/ba3f.py)
- [x] [BA3G: Find an Eulerian Path in a Graph](rosalind/bioinformatics_textbook_track/ba3g.py)
- [x] [BA3H: Reconstruct a String from its k-mer Composition](rosalind/bioinformatics_textbook_track/ba3h.py)
- [x] [BA3I: Find a k-Universal Circular String](rosalind/bioinformatics_textbook_track/ba3i.py)
- [x] [BA3J: Reconstruct a String from its Paired Composition](rosalind/bioinformatics_textbook_track/ba3j.py)
- [x] [BA3K: Generate Contigs from a Collection of Reads](rosalind/bioinformatics_textbook_track/ba3k.py)
- [x] [BA3L: Construct a String Spelled by a Gapped Genome Path](rosalind/bioinformatics_textbook_track/ba3l.py)
- [x] [BA3M: Generate All Maximal Non-Branching Paths in a Graph](rosalind/bioinformatics_textbook_track/ba3m.py)
- [x] [BA4A: Translate an RNA String into an Amino Acid String](rosalind/bioinformatics_textbook_track/ba4a.py)
- [x] [BA4B: Find Substrings of a Genome Encoding a Given Amino Acid String](rosalind/bioinformatics_textbook_track/ba4b.py)
- [x] [BA4C: Generate the Theoretical Spectrum of a Cyclic Peptide](rosalind/bioinformatics_textbook_track/ba4c.py)
- [x] [BA4D: Compute the Number of Peptides of Given Total Mass](rosalind/bioinformatics_textbook_track/ba4d.py)
- [x] [BA4E: Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum](rosalind/bioinformatics_textbook_track/ba4e.py)
- [x] [BA4F: Compute the Score of a Cyclic Peptide Against a Spectrum](rosalind/bioinformatics_textbook_track/ba4f.py)
- [x] [BA4G: Implement LeaderboardCyclopeptideSequencing](rosalind/bioinformatics_textbook_track/ba4g.py)
- [x] [BA4H: Generate the Convolution of a Spectrum](rosalind/bioinformatics_textbook_track/ba4h.py)
- [x] [BA4I: Implement ConvolutionCyclopeptideSequencing](rosalind/bioinformatics_textbook_track/ba4i.py)
- [x] [BA4J: Generate the Theoretical Spectrum of a Linear Peptide](rosalind/bioinformatics_textbook_track/ba4j.py)
- [x] [BA4K: Compute the Score of a Linear Peptide](rosalind/bioinformatics_textbook_track/ba4k.py)
- [x] [BA4L: Trim a Peptide Leaderboard](rosalind/bioinformatics_textbook_track/ba4l.py)
- [x] [BA4M: Solve the Turnpike Problem](rosalind/bioinformatics_textbook_track/ba4m.py)
- [x] [BA5A: Find the Minimum Number of Coins Needed to Make Change](rosalind/bioinformatics_textbook_track/ba5a.py)
- [x] [BA5B: Find the Length of a Longest Path in a Manhattan-like Grid](rosalind/bioinformatics_textbook_track/ba5b.py)
- [x] [BA5C: Find a Longest Common Subsequence of Two Strings](rosalind/bioinformatics_textbook_track/ba5c.py)
- [x] [BA5D: Find the Longest Path in a DAG](rosalind/bioinformatics_textbook_track/ba5d.py)
- [x] [BA5E: Find a Highest-Scoring Alignment of Two Strings](rosalind/bioinformatics_textbook_track/ba5e.py)
- [x] [BA5F: Find a Highest-Scoring Local Alignment of Two Strings](rosalind/bioinformatics_textbook_track/ba5f.py)
- [x] [BA5G: Compute the Edit Distance Between Two Strings](rosalind/bioinformatics_textbook_track/ba5g.py)
- [x] [BA5H: Find a Highest-Scoring Fitting Alignment of Two Strings](rosalind/bioinformatics_textbook_track/ba5h.py)
- [x] [BA5I: Find a Highest-Scoring Overlap Alignment of Two Strings](rosalind/bioinformatics_textbook_track/ba5i.py)
- [x] [BA5J: Align Two Strings Using Affine Gap Penalties](rosalind/bioinformatics_textbook_track/ba5j.py)
- [x] [BA5K: Find a Middle Edge in an Alignment Graph in Linear Space](rosalind/bioinformatics_textbook_track/ba5k.py)
- [x] [BA5L: Align Two Strings Using Linear Space](rosalind/bioinformatics_textbook_track/ba5l.py)
- [x] [BA5M: Find a Highest-Scoring Multiple Sequence Alignment](rosalind/bioinformatics_textbook_track/ba5m.py)
- [x] [BA5N: Find a Topological Ordering of a DAG](rosalind/bioinformatics_textbook_track/ba5n.py)
- [x] [BA6A: Implement GreedySorting to Sort a Permutation by Reversals](rosalind/bioinformatics_textbook_track/ba6a.py)
- [x] [BA6B: Compute the Number of Breakpoints in a Permutation](rosalind/bioinformatics_textbook_track/ba6b.py)
- [x] [BA6C: Compute the 2-Break Distance Between a Pair of Genomes](rosalind/bioinformatics_textbook_track/ba6c.py)
- [x] [BA6D: Find a Shortest Transformation of One Genome into Another by 2-Breaks](rosalind/bioinformatics_textbook_track/ba6d.py)
- [x] [BA6E: Find All Shared k-mers of a Pair of Strings](rosalind/bioinformatics_textbook_track/ba6e.py)
- [x] [BA6F: Implement ChromosomeToCycle](rosalind/bioinformatics_textbook_track/ba6f.py)
- [x] [BA6G: Implement CycleToChromosome](rosalind/bioinformatics_textbook_track/ba6g.py)
- [x] [BA6H: Implement ColoredEdges](rosalind/bioinformatics_textbook_track/ba6h.py)
- [x] [BA6I: Implement GraphToGenome](rosalind/bioinformatics_textbook_track/ba6i.py)
- [x] [BA6J: Implement 2-BreakOnGenomeGraph](rosalind/bioinformatics_textbook_track/ba6j.py)
- [x] [BA6K: Implement 2-BreakOnGenome](rosalind/bioinformatics_textbook_track/ba6k.py)
- [x] [BA7A: Compute Distances Between Leaves](rosalind/bioinformatics_textbook_track/ba7a.py)
- [x] [BA7B: Compute Limb Lengths in a Tree](rosalind/bioinformatics_textbook_track/ba7b.py)
- [x] [BA7C: Implement AdditivePhylogeny](rosalind/bioinformatics_textbook_track/ba7c.py)
- [x] [BA7D: Implement UPGMA](rosalind/bioinformatics_textbook_track/ba7d.py)
- [x] [BA7E: Implement the Neighbor Joining Algorithm](rosalind/bioinformatics_textbook_track/ba7e.py)
- [x] [BA7F: Implement SmallParsimony](rosalind/bioinformatics_textbook_track/ba7f.py)
- [x] [BA7G: Adapt SmallParsimony to Unrooted Trees](rosalind/bioinformatics_textbook_track/ba7g.py)
- [x] [BA8A: Implement FarthestFirstTraversal](rosalind/bioinformatics_textbook_track/ba8a.py)
- [x] [BA8B: Compute the Squared Error Distortion](rosalind/bioinformatics_textbook_track/ba8ab.py)
- [x] [BA8C: Implement the Lloyd Algorithm for k-Means Clustering](rosalind/bioinformatics_textbook_track/ba8c.py)
- [x] [BA8D: Implement the Soft k-Means Clustering Algorithm](rosalind/bioinformatics_textbook_track/ba8d.py)
- [x] [BA8E: Implement Hierarchical Clustering](rosalind/bioinformatics_textbook_track/ba8e.py)
- [x] [BA9A: Construct a Trie from a Collection of Patterns](rosalind/bioinformatics_textbook_track/ba9a.py)
- [x] [BA9B: Implement TrieMatching](rosalind/bioinformatics_textbook_track/ba9b.py)
- [x] [BA9C: Construct the Suffix Tree of a String](rosalind/bioinformatics_textbook_track/ba9c.py)
- [x] [BA9D: Find the Longest Repeat in a String](rosalind/bioinformatics_textbook_track/ba9d.py)
- [x] [BA9E: Find the Longest Substring Shared by Two Strings](rosalind/bioinformatics_textbook_track/ba9e.py)
- [x] [BA9F: Find the Shortest Non-Shared Substring of Two Strings](rosalind/bioinformatics_textbook_track/ba9f.py)
- [x] [BA9G: Construct the Suffix Array of a String](rosalind/bioinformatics_textbook_track/ba9g.py)
- [x] [BA9H: Pattern Matching with the Suffix Array](rosalind/bioinformatics_textbook_track/ba9h.py)
- [x] [BA9I: Construct the Burrows-Wheeler Transform of a String](rosalind/bioinformatics_textbook_track/ba9i.py)
- [x] [BA9J: Reconstruct a String from its Burrows-Wheeler Transform](rosalind/bioinformatics_textbook_track/ba9j.py)
- [x] [BA9K: Generate the Last-to-First Mapping of a String](rosalind/bioinformatics_textbook_track/ba9k.py)
- [x] [BA9L: Implement BWMatching](rosalind/bioinformatics_textbook_track/ba9l.py)
- [x] [BA9M: Implement BetterBWMatching](rosalind/bioinformatics_textbook_track/ba9m.py)
- [ ] BA9N: Find All Occurrences of a Collection of Patterns in a String
- [ ] BA9O: Find All Approximate Occurrences of a Collection of Patterns in a String
- [x] [BA9P: Implement TreeColoring](rosalind/bioinformatics_textbook_track/ba9p.py)
- [x] [BA9Q: Construct the Partial Suffix Array of a String](rosalind/bioinformatics_textbook_track/ba9q.py)
- [ ] BA9R: Construct a Suffix Tree from a Suffix Array
- [x] [BA10A: Compute the Probability of a Hidden Path](rosalind/bioinformatics_textbook_track/ba10a.py)
- [x] [BA10B: Compute the Probability of an Outcome Given a Hidden Path](rosalind/bioinformatics_textbook_track/ba10b.py)
- [x] [BA10C: Implement the Viterbi Algorithm](rosalind/bioinformatics_textbook_track/ba10c.py)
- [x] [BA10D: Compute the Probability of a String Emitted by an HMM](rosalind/bioinformatics_textbook_track/ba10d.py)
- [x] [BA10E: Construct a Profile HMM](rosalind/bioinformatics_textbook_track/ba10e.py)
- [x] [BA10F: Construct a Profile HMM with Pseudocounts](rosalind/bioinformatics_textbook_track/ba10f.py)
- [ ] BA10G: Perform a Multiple Sequence Alignment with a Profile HMM
- [x] [BA10H: Estimate the Parameters of an HMM](rosalind/bioinformatics_textbook_track/ba10h.py)
- [x] [BA10I: Implement Viterbi Learning](rosalind/bioinformatics_textbook_track/ba10i.py)
- [x] [BA10J: Solve the Soft Decoding Problem](rosalind/bioinformatics_textbook_track/ba10j.py)
- [x] [BA10K: Implement Baum-Welch Learning](rosalind/bioinformatics_textbook_track/ba10j.py)
- [x] [BA11A: Construct the Graph of a Spectrum](rosalind/bioinformatics_textbook_track/ba11a.py)
- [x] [BA11B: Implement DecodingIdealSpectrum](rosalind/bioinformatics_textbook_track/ba11b.py)
- [x] [BA11C: Convert a Peptide into a Peptide Vector](rosalind/bioinformatics_textbook_track/ba11c.py)
- [x] [BA11D: Convert a Peptide Vector into a Peptide](rosalind/bioinformatics_textbook_track/ba11d.py)
- [x] [BA11E: Sequence a Peptide](rosalind/bioinformatics_textbook_track/ba11e.py)
- [x] [BA11F: Find a Highest-Scoring Peptide in a Proteome against a Spectrum](rosalind/bioinformatics_textbook_track/ba11f.py)
- [x] [BA11G: Implement PSMSearch](rosalind/bioinformatics_textbook_track/ba11g.py)
- [x] [BA11H: Compute the Size of a Spectral Dictionary](rosalind/bioinformatics_textbook_track/ba11h.py)
- [x] [BA11I: Compute the Probability of a Spectral Dictionary](rosalind/bioinformatics_textbook_track/ba11i.py)
- [ ] BA11J: Find a Highest-Scoring Modified Peptide against a Spectrum

### Algorithmic Heights

- [x] [FIBO: Fibonacci Numbers](rosalind/algorithmic_heights/fibo.py)
- [x] [BINS: Binary Search](rosalind/algorithmic_heights/bins.py)
- [x] [DEG: Degree Array](rosalind/algorithmic_heights/deg.py)
- [x] [INS: Insertion Sort](rosalind/algorithmic_heights/ins.py)
- [x] [DDEG: Double-Degree Array](rosalind/algorithmic_heights/ddeg.py)
- [x] [MAJ: Majority Element](rosalind/algorithmic_heights/maj.py)
- [x] [MER: Merge Two Sorted Arrays](rosalind/algorithmic_heights/mer.py)
- [x] [2SUM: 2SUM](rosalind/algorithmic_heights/2sum.py)
- [x] [BFS: Breadth-First Search](rosalind/algorithmic_heights/bfs.py)
- [x] [CC: Connected Components](rosalind/algorithmic_heights/cc.py)
- [x] [MS: Merge Sort](rosalind/algorithmic_heights/ms.py)
- [x] [PAR: 2-Way Partition](rosalind/algorithmic_heights/par.py)
- [x] [3SUM: 3SUM](rosalind/algorithmic_heights/3sum.py)
- [x] [DAG: Testing Acyclicity](rosalind/algorithmic_heights/dag.py)
- [x] [DIJ: Dijkstra's Algorithm](rosalind/algorithmic_heights/dij.py)
- [x] [INV: Counting Inversions](rosalind/algorithmic_heights/inv.py)
- [x] [TS: Topological Sorting](rosalind/algorithmic_heights/ts.py)
- [x] [BIP: Testing Bipartiteness](rosalind/algorithmic_heights/bip.py)
- [x] [PAR3: 3-Way Partition](rosalind/algorithmic_heights/par3.py)
- [x] [MED: Median](rosalind/algorithmic_heights/med.py)
- [x] [QS: Quick Sort](rosalind/algorithmic_heights/qs.py)
- [x] [HEA: Building a Heap](rosalind/algorithmic_heights/hea.py)
- [ ] SQ: Square in a Graph
- [x] [BF: Bellman-Ford Algorithm](rosalind/algorithmic_heights/bf.py)
- [x] [CTE: Shortest Cycle Through a Given Edge](rosalind/algorithmic_heights/cte.py)
- [x] [HS: Heap Sort](rosalind/algorithmic_heights/hs.py)
- [x] [PS: Partial Sort](rosalind/algorithmic_heights/ps.py)
- [x] [HDAG: Hamiltonian Path in DAG](rosalind/algorithmic_heights/hdag.py)
- [x] [NWC: Negative Weight Cycle](rosalind/algorithmic_heights/nwc.py)
- [x] [SCC: Strongly Connected Components](rosalind/algorithmic_heights/scc.py)
- [ ] 2SAT: 2-Satisfiability
- [ ] GS: General Sink
- [ ] SC: Semi-Connected Graph
- [ ] SDAG: Shortest Paths in DAG


[rosalind.info]: https://rosalind.info
[locations]: https://rosalind.info/problems/locations/
[Python Village]: https://rosalind.info/problems/list-view/?location=python_village
[Bioinformatics Stronghold]: https://rosalind.info/problems/list-view/
[Algorithmic Heights]: https://rosalind.info/problems/list-view/?location=algorithmic_heights
[Bioinformatics Armory]: https://rosalind.info/problems/list-view/?location=bioinformatics_armory
[Bioinformatics Textbook Track]: https://rosalind.info/problems/list-view/?location=bioinformatics_textbook_track
[Bioinformatics Algorithms: An Active Learning Approach]: https://www.bioinformaticsalgorithms.org/
[Algorithms]: https://www.google.co.uk/books/edition/Algorithms/DJSUCgAAQBAJ?hl=en
[poetry]: https://python-poetry.org/
[typer]: https://typer.tiangolo.com/
