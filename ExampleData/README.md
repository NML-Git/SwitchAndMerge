# Example data

This folder contains small example input and output files

# Example run of PairedEnd_SwitchAndMerge.m

The two fastw files `Testset_fastq_1` and `Testset_fastq_2` contains a small test set of 10,000 reads, with mixed forward and reverse reads.

When run through the `PairedEnd_SwitchAndMerge.m` script they will generate output file corresponding to `Output_TestSet_SwitchedAndMerged.fastq` and `Output_TestSet_SwitchedAndMerged_log.txt`.

The user may compare their log-file with the example log `Output_TestSet_SwitchedAndMerged_log.txt`
	Paths, filenames and timestamps have been anonymized.

The script was run with default parameters, except for the output-specific parameter `Show_every_nth_read` which was set to 1,000 to accommodate the small test-set.
