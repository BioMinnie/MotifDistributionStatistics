Motif Distribution Statistics
=============================

The following scripts are tailored to suite the data published in Forde et al. 2015
Running title: *E. coli* EC958 methylome

Python script dedicated to counting motifs within genomic segments
R scripts dedicated to the statistical analysis and graphical output of the distribution of methylation motifs

Scripts available:

  # motif_counts.py
  
  # adjusted_anova.R
  
  # motif_distance_vs_position.R

How to:

motif_analysis.py

  1. Takes in a tab delimited file of motif positions in genome
  2. Takes in a tab delimited file of segments (in this case, 1 kbp segments with a 250bp overlap of the genome, generated using Bedtools makewindows option)
  3. Counts number of motifs within each segment and writes output to the outfile

adjusted_anova.R

  1. Takes in a tab delimited file, with each motif count per segment and region of each segment required
  2. For each methylation motif in the genome, perform an analysis of variance, adjusted for heteroscedasticity. Reports anova      test statistic, mean motif count per genomic region, pair-wise post hoc P values and confidence intervals around the mean

motif_distance_vs_position.R

  1. Takes in a tab delimited file (start position of motif, distance between motifs, genomic region (i.e. Core-genome, Prophage, Genomic Island)), computes summary statistics and plots the distance between motifs (bp) vs. position in genome (bp). Figure can be saved in svg format
