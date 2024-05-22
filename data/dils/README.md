This repository contains the data files used for demographic inference (DILS) analysis, details about the input format for the DILS software can be found at the following link [https://github.com/dils-popgen/dils/blob/master/manual.pdf](https://github.com/dils-popgen/dils/blob/master/manual.pdf)

The *.yaml* file is a configuration file that stores data paths, model & population parameters for the DILS analysis 
The fasta alignment file combines 935 random non-coding hamlet genome sequences, separated in 2 alleles for each of the 335 genomes (20 hamlet species + 3 outgroup species).

Content tree:
<pre>
.
├── large_small_ser_2pop.yaml # configuration file input for DILS
└── all935.fa  # fasta alignment file
</pre>
