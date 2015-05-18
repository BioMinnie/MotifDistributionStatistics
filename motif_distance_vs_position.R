#######################################################################################################
#                                                                                                     #
## AUTHOR: Melinda M Ashcroft                                                                         #
## AFFILIATION: Beatson Lab | SCMB - University of Queensland                                         #
#                                                                                                     #
## PURPOSE: Computes summary statistics and plots the distance between motifs (bp) vs. position       #
## Taking into account core-genome, Prophages and Genomic Islands                                     #
#                                                                                                     #
## PUBLISHED: The following script is tailored to suit the data published in Forde et al. 2015        #
#######################################################################################################

# Import library
library(ggplot2)
library(scales)
library(RSvgDevice)

# Read in the data
motif_reg = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(motif_reg)

# Check first few observations
head(motif_reg)

# Summary statistics
dist_sd <- sd(motif_reg$distance)
dist_mean <- mean(motif_reg$distance)

# Generate cutoff
cutoff <- dist_mean + (3 * dist_sd)
cutoff
hist(motif_reg$distance)

# Normal scatterplot - no colour
distance = motif_reg$distance
position = motif_reg$start
plot(position, distance)

# Scatterplot using own colours
col_plot <- ggplot(motif_reg, aes(start, distance, color=factor(region))) + 
  xlab("Position in EC958 (bp)") +
  ylab("Distance between GATC motifs (bp)") +
  geom_point() + scale_color_manual(values=c("#BABABA", "#4169E1", "#DB3D95"), 
                                    name = "Region",
                                    breaks=c("genome", "Genomic_Island", "Prophage"), 
                                    labels=c("Genome", "Genomic Island", "Prophage")) + 
  geom_hline(yintercept=1119.55, linetype="dashed") +
  theme_bw(base_size = 12, base_family = "Helvetica") + theme(panel.border = element_rect(color = "dark grey", size = .3), 
                                                             panel.grid.major.y = element_line(size = .2, color = "grey"), 
                                                             panel.grid.minor.y = element_line(size = .2, color = "#E8E7E7"),
                                                             panel.grid.major.x = element_line(size = .2, color = "#E8E7E7"), 
                                                             panel.grid.minor.x = element_blank(),
                                                             axis.line = element_line(size=.7, color = "black"),
                                                             text = element_text(size=14)) 
# Generate plot
col_plot

# Save the plot as an svg image file
devSVG(file="motifs_distribution.svg", width=7, height=3.85) # Open plotting device
col_plot
dev.off() # Turn off plotting device
