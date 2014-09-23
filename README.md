Motif Distribution Statistics
=============================

Python script dedicated to counting motifs within DNA sequences.

R scripts dedicated to the statistical analysis of the distribution of methylation motifs.

Scripts available:

  motif_analysis.py
  adjusted_anova.R
  motif_distance_vs_position.R

How to:

motif_analysis.py

1. Takes in a tab delimited file of motif positions in genome
2. Takes in a tab delimited file of 1kb segments with a 250bp overlap of the genome, generated using bedtools makewindows option
3. Counts number of motifs within each 1kb segment and reports counts, segment start and end positions, segment name and genomic region (genome, mobile genetic elements) and writes output to the outfile

adjusted_anova.R

1. Takes in a tab delimited file (each motif count per segment and region required)
2. For each methylation motif in the genome, perform an analysis of variance, adjusted for heteroscedasticity. Reports anova test statistic, mean motif count per genomic region, pair-wise post hoc P values and confidence intervals around the mean.

motif_distance_vs_position.R

1. Takes in a tab delimited file (start position of motif, distance between motifs, genomic region), computes summary statistics and plots the distance between motifs (bp) vs. position in genome (bp). Figure can be saved in either svg or pdf format.
