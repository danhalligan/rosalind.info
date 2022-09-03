# Bioinformatics Armory Solutions

Solutions to "[Bioinformatics Armory]" coding challenges from [rosalind.info].

The solutions here are python based only. For some solutions you'll need
to use installed bioinformatics tools (or web interfaces).

[Bioinformatics Armory]: https://rosalind.info/problems/list-view/?location=bioinformatics-armory
[rosalind.info]: https://rosalind.info

## Running the solutions

Each problem is solved in a separate file.

To run this tool to solve a particular problem, you can either run the script
(providing the input file as a command line argument) or use `solve.py` which
will solve based on the input file name.

```{shell}
python3 ini.py test/rosalind_ini.txt
```

```{shell}
./solve.py test/rosalind_ini.txt
```

## Notes

* For "meme" and "clus" there are solutions. Use the web interface as
  instructed.
* For "subo", you need to run the online interface, identify the 32-40 bp and
  then can use the solution provided to count the occurrences of this in
  the sequences.
