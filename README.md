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
poetry run rosalind ini2 ~/Downloads/rosalind_ini2.txt
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

- [x] [INI1: Installing Python](python-village/ini1.py)
- [x] [INI2: Variables and Some Arithmetic](python-village/ini2.py)
- [x] [INI3: Strings and Lists](python-village/ini3.py)
- [x] [INI4: Conditions and Loops](python-village/ini4.py)
- [x] [INI5: Working with Files](python-village/ini5.py)
- [x] [INI6: Dictionaries](python-village/ini6.py)

### Bioinformatics Armory

- [x] [INI: Introduction to the Bioinformatics Armory](ini.py)
- [x] [GBK: GenBank Introduction](gbk.py)
- [x] [FRMT: Data Formats](frmt.py)
- [x] MEME: New Motif Discovery
- [x] [NEED: Pairwise Global Alignment](need.py)
- [x] [TFSQ: FASTQ format introduction](tfsq.py)
- [x] [PHRE: Read Quality Distribution](phre.py)
- [x] [PTRA: Protein Translation](ptra.py)
- [x] [FILT: Read Filtration by Quality](filt.py)
- [x] [RVCO: Complementing a Strand of DNA](rvco.py)
- [x] [SUBO: Suboptimal Local Alignment](subo.py)
- [x] [BPHR: Base Quality Distribution](bphr.py)
- [x] CLUS: Global Multiple Alignment
- [x] [ORFR: Finding Genes with ORFs](orfr.py)
- [x] [BFIL: Base Filtration by Quality](bfil.py)

#### Notes

* For "MEME" and "CLUS" I have not written a solution. Use the web interface as
  instructed.
* For "SUBO", you need to run the online interface, identify the 32-40 bp and
  then can use the solution here to count the occurrences of this in
  the sequences.

### Bioinformatics Textbook Track

- [x] [BA1A: Compute the Number of Times a Pattern Appears in a Text](ba1a.py)
- [x] [BA1B: Find the Most Frequent Words in a String](ba1b.py)
- [x] [BA1C: Find the Reverse Complement of a String](ba1c.py)
- [x] [BA1D: Find All Occurrences of a Pattern in a String](ba1d.py)
- [x] [BA1E: Find Patterns Forming Clumps in a String](ba1e.py)
- [x] [BA1F: Find a Position in a Genome Minimizing the Skew](ba1f.py)
- [x] [BA1G: Compute the Hamming Distance Between Two Strings](ba1g.py)
- [x] [BA1H: Find All Approximate Occurrences of a Pattern in a String](ba1h.py)
- [x] [BA1I: Find the Most Frequent Words with Mismatches in a String](ba1i.py)
- [x] [BA1J: Find Frequent Words with Mismatches and Reverse Complements](ba1j.py)
- [x] [BA1K: Generate the Frequency Array of a String](ba1k.py)
- [x] [BA1L: Implement PatternToNumber](ba1l.py)
- [x] [BA1M: Implement NumberToPattern](ba1m.py)
- [x] [BA1N: Generate the d-Neighborhood of a String](ba1n.py)
- [x] [BA2A: Implement MotifEnumeration](ba2a.py)
- [x] [BA2B: Find a Median String](ba2b.py)
- [x] [BA2C: Find a Profile-most Probable k-mer in a String](ba2c.py)
- [x] [BA2D: Implement GreedyMotifSearch](ba2d.py)
- [x] [BA2E: Implement GreedyMotifSearch with Pseudocounts](ba2e.py)
- [x] [BA2F: Implement RandomizedMotifSearch](ba2f.py)
- [x] [BA2G: Implement GibbsSampler](ba2g.py)
- [x] [BA2H: Implement DistanceBetweenPatternAndStrings](ba2h.py)
- [x] [BA3A: Generate the k-mer Composition of a String](ba3a.py)
- [x] [BA3B: Reconstruct a String from its Genome Path](ba3b.py)
- [x] [BA3C: Construct the Overlap Graph of a Collection of k-mers](ba3c.py)
- [x] [BA3D: Construct the De Bruijn Graph of a String](ba3d.py)
- [x] [BA3E: Construct the De Bruijn Graph of a Collection of k-mers](ba3e.py)
- [x] [BA3F: Find an Eulerian Cycle in a Graph](ba3f.py)
- [x] [BA3G: Find an Eulerian Path in a Graph](ba3g.py)
- [x] [BA3H: Reconstruct a String from its k-mer Composition](ba3h.py)
- [x] [BA3I: Find a k-Universal Circular String](ba3i.py)
- [x] [BA3J: Reconstruct a String from its Paired Composition](ba3j.py)
- [x] [BA3K: Generate Contigs from a Collection of Reads](ba3k.py)
- [x] [BA3L: Construct a String Spelled by a Gapped Genome Path](ba3l.py)
- [x] [BA3M: Generate All Maximal Non-Branching Paths in a Graph](ba3m.py)
- [x] [BA4A: Translate an RNA String into an Amino Acid String](ba4a.py)
- [x] [BA4B: Find Substrings of a Genome Encoding a Given Amino Acid String](ba4b.py)
- [x] [BA4C: Generate the Theoretical Spectrum of a Cyclic Peptide](ba4c.py)
- [x] [BA4D: Compute the Number of Peptides of Given Total Mass](ba4d.py)
- [x] [BA4E: Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum](ba4e.py)
- [x] [BA4F: Compute the Score of a Cyclic Peptide Against a Spectrum](ba4f.py)
- [x] [BA4G: Implement LeaderboardCyclopeptideSequencing](ba4g.py)
- [x] [BA4H: Generate the Convolution of a Spectrum](ba4h.py)
- [x] [BA4I: Implement ConvolutionCyclopeptideSequencing](ba4i.py)
- [x] [BA4J: Generate the Theoretical Spectrum of a Linear Peptide](ba4j.py)
- [x] [BA4K: Compute the Score of a Linear Peptide](ba4k.py)
- [x] [BA4L: Trim a Peptide Leaderboard](ba4l.py)
- [x] [BA4M: Solve the Turnpike Problem](ba4m.py)
- [x] [BA5A: Find the Minimum Number of Coins Needed to Make Change](ba5a.py)
- [x] [BA5B: Find the Length of a Longest Path in a Manhattan-like Grid](ba5b.py)
- [x] [BA5C: Find a Longest Common Subsequence of Two Strings](ba5c.py)
- [x] [BA5D: Find the Longest Path in a DAG](ba5d.py)
- [x] [BA5E: Find a Highest-Scoring Alignment of Two Strings](ba5e.py)
- [x] [BA5F: Find a Highest-Scoring Local Alignment of Two Strings](ba5f.py)
- [x] [BA5G: Compute the Edit Distance Between Two Strings](ba5g.py)
- [x] [BA5H: Find a Highest-Scoring Fitting Alignment of Two Strings](ba5h.py)
- [x] [BA5I: Find a Highest-Scoring Overlap Alignment of Two Strings](ba5i.py)
- [x] [BA5J: Align Two Strings Using Affine Gap Penalties](ba5j.py)
- [x] [BA5K: Find a Middle Edge in an Alignment Graph in Linear Space](ba5k.py)
- [x] [BA5L: Align Two Strings Using Linear Space](ba5l.py)
- [x] [BA5M: Find a Highest-Scoring Multiple Sequence Alignment](ba5m.py)
- [x] [BA5N: Find a Topological Ordering of a DAG](ba5n.py)
- [x] [BA6A: Implement GreedySorting to Sort a Permutation by Reversals](ba6a.py)
- [x] [BA6B: Compute the Number of Breakpoints in a Permutation](ba6b.py)
- [x] [BA6C: Compute the 2-Break Distance Between a Pair of Genomes](ba6c.py)
- [x] [BA6D: Find a Shortest Transformation of One Genome into Another by 2-Breaks](ba6d.py)
- [x] [BA6E: Find All Shared k-mers of a Pair of Strings](ba6e.py)
- [x] [BA6F: Implement ChromosomeToCycle](ba6f.py)
- [x] [BA6G: Implement CycleToChromosome](ba6g.py)
- [x] [BA6H: Implement ColoredEdges](ba6h.py)
- [x] [BA6I: Implement GraphToGenome](ba6i.py)
- [x] [BA6J: Implement 2-BreakOnGenomeGraph](ba6j.py)
- [x] [BA6K: Implement 2-BreakOnGenome](ba6k.py)
- [x] [BA7A: Compute Distances Between Leaves](ba7a.py)
- [x] [BA7B: Compute Limb Lengths in a Tree](ba7b.py)
- [x] [BA7C: Implement AdditivePhylogeny](ba7c.py)
- [ ] BA7D: Implement UPGMA
- [ ] BA7E: Implement the Neighbor Joining Algorithm
- [ ] BA7F: Implement SmallParsimony
- [ ] BA7G: Adapt SmallParsimony to Unrooted Trees
- [ ] BA8A: Implement FarthestFirstTraversal
- [ ] BA8B: Compute the Squared Error Distortion
- [ ] BA8C: Implement the Lloyd Algorithm for k-Means Clustering
- [ ] BA8D: Implement the Soft k-Means Clustering Algorithm
- [ ] BA8E: Implement Hierarchical Clustering
- [ ] BA9A: Construct a Trie from a Collection of Patterns
- [ ] BA9B: Implement TrieMatching
- [ ] BA9C: Construct the Suffix Tree of a String
- [ ] BA9D: Find the Longest Repeat in a String
- [ ] BA9E: Find the Longest Substring Shared by Two Strings
- [ ] BA9F: Find the Shortest Non-Shared Substring of Two Strings
- [ ] BA9G: Construct the Suffix Array of a String
- [ ] BA9H: Pattern Matching with the Suffix Array
- [ ] BA9I: Construct the Burrows-Wheeler Transform of a String
- [ ] BA9J: Reconstruct a String from its Burrows-Wheeler Transform
- [ ] BA9K: Generate the Last-to-First Mapping of a String
- [ ] BA9L: Implement BWMatching
- [ ] BA9M: Implement BetterBWMatching
- [ ] BA9N: Find All Occurrences of a Collection of Patterns in a String
- [ ] BA9O: Find All Approximate Occurrences of a Collection of Patterns in a String
- [ ] BA9P: Implement TreeColoring
- [ ] BA9Q: Construct the Partial Suffix Array of a String
- [ ] BA9R: Construct a Suffix Tree from a Suffix Array
- [ ] BA10A: Compute the Probability of a Hidden Path
- [ ] BA10B: Compute the Probability of an Outcome Given a Hidden Path
- [ ] BA10C: Implement the Viterbi Algorithm
- [ ] BA10D: Compute the Probability of a String Emitted by an HMM
- [ ] BA10E: Construct a Profile HMM
- [ ] BA10F: Construct a Profile HMM with Pseudocounts
- [ ] BA10G: Perform a Multiple Sequence Alignment with a Profile HMM
- [ ] BA10H: Estimate the Parameters of an HMM
- [ ] BA10I: Implement Viterbi Learning
- [ ] BA10J: Solve the Soft Decoding Problem
- [ ] BA10K: Implement Baum-Welch Learning
- [ ] BA11A: Construct the Graph of a Spectrum
- [ ] BA11B: Implement DecodingIdealSpectrum
- [ ] BA11C: Convert a Peptide into a Peptide Vector
- [ ] BA11D: Convert a Peptide Vector into a Peptide
- [ ] BA11E: Sequence a Peptide
- [ ] BA11F: Find a Highest-Scoring Peptide in a Proteome against a Spectrum
- [ ] BA11G: Implement PSMSearch
- [ ] BA11H: Compute the Size of a Spectral Dictionary
- [ ] BA11I: Compute the Probability of a Spectral Dictionary
- [ ] BA11J: Find a Highest-Scoring Modified Peptide against a Spectrum

### Algorithmic Heights

- [x] [FIBO: Fibonacci Numbers](fibo.py)
- [x] [BINS: Binary Search](bins.py)
- [x] [DEG: Degree Array](deg.py)
- [x] [INS: Insertion Sort](ins.py)
- [x] [DDEG: Double-Degree Array](ddeg.py)
- [x] [MAJ: Majority Element](maj.py)
- [x] [MER: Merge Two Sorted Arrays](mer.py)
- [x] [2SUM: 2SUM](2sum.py)
- [x] [BFS: Breadth-First Search](bfs.py)
- [x] [CC: Connected Components](cc.py)
- [x] [MS: Merge Sort](ms.py)
- [x] [PAR: 2-Way Partition](par.py)
- [x] [3SUM: 3SUM](3sum.py)
- [x] [DAG: Testing Acyclicity](dag.py)
- [x] [DIJ: Dijkstra's Algorithm](dij.py)
- [x] [INV: Counting Inversions](inv.py)
- [x] [TS: Topological Sorting](ts.py)
- [ ] HEA: Building a Heap
- [ ] BIP: Testing Bipartiteness
- [ ] PAR3: 3-Way Partition
- [ ] SQ: Square in a Graph
- [ ] BF: Bellman-Ford Algorithm
- [ ] CTE: Shortest Cycle Through a Given Edge
- [ ] HS: Heap Sort
- [ ] MED: Median
- [ ] PS: Partial Sort
- [ ] HDAG: Hamiltonian Path in DAG
- [ ] NWC: Negative Weight Cycle
- [ ] QS: Quick Sort
- [ ] SCC: Strongly Connected Components
- [ ] 2SAT: 2-Satisfiability
- [ ] GS: General Sink
- [ ] SC: Semi-Connected Graph
- [ ] SDAG: Shortest Paths in DAG


[rosalind.info]: https://rosalind.info
[locations]: https://rosalind.info/problems/locations/
[Python Village]: https://rosalind.info/problems/list-view/?location=python-village
[Bioinformatics Stronghold]: https://rosalind.info/problems/list-view/
[Algorithmic Heights]: https://rosalind.info/problems/list-view/?location=algorithmic-heights
[Bioinformatics Armory]: https://rosalind.info/problems/list-view/?location=bioinformatics-armory
[Bioinformatics Textbook Track]: https://rosalind.info/problems/list-view/?location=bioinformatics-textbook-track
[Bioinformatics Algorithms: An Active Learning Approach]: https://www.bioinformaticsalgorithms.org/
[Algorithms]: https://www.google.co.uk/books/edition/Algorithms/DJSUCgAAQBAJ?hl=en
[poetry]: https://python-poetry.org/
[typer]: https://typer.tiangolo.com/
