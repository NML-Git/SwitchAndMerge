# Post-processing scripts

These scripts perform downstream analysis of library data that has been subjected to PHASTpep pt 1, to get an output of unique sequences and counts of each. Some of the scripts may also take as input the PHASTpep pt 2 output. The exact formats needed are specified in the header of the individual scripts and the sections below.

The scripts generate a workflow that allows one to remove known target-unrelated peptides (TUPs) from a given library with `RemoveTUPsFromCounts.m`, after which we can locate sequence motif patterns that still remains in the library with `SearchForRepeats.m`. The exact matches to specific motifs may be elucidated with `ShowRepeatResults.m`.

Finally, the diversity profile of a library (before or after normalization as described in the folder "PHASTpep compatibility") may be illustrated using `PlotDiversity.m`. Diversity profiles are a useful tool to inspect the large sequence library dataset.

| Script | Purpose |
|----------|----------|
| `RemoveTUPsFromCounts.m` | Completely removes sequences that partially match a given set of sequences |
| `SearchForRepeats.m` | Finds motifs of sequence families within a library |
| `ShowRepeatResults.m` | Outputs sequences in a library that match a given set of motifs |
| `PlotDiversity.m` | Plots diversity profile of libraries |

# Typical workflow
1. Remove known TUPs using `RemoveTUPsFromCounts.m`.
2. Identify recurring sequence motifs with `SearchForRepeats.m`.
3. Inspect motif-containing sequences with `ShowRepeatResults.m`.
4. Visualize library diversity using `PlotDiversity.m`.

Detailed input/output specifications and links to example data are provided in the header of each script.
