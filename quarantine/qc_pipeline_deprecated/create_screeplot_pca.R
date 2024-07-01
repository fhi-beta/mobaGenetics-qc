#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

### INPUT ARGS ###

# Input is the .pca output from Eigenstrat
eigval.input <- args[1]
#eigval.input <- '/home/oyvind/harvest/media/local-disk/oyvindhe/erc-genotypes-results/moba24-erc-genotypes/offspring-qc/data/pca-final/final_data_core_samples_only.pca'

# number of PCAs in file
num.pca <- as.numeric(args[2])
  
# Header line for plot
plot.header <- args[3]

# output path of plot
out.path <- args[4]

# output name of plot pdf
out.plotname <- args[5]

### SCRIPT ###
library(ggplot2)

# load PCA file
pca <- read.table(eigval.input, header=F, nrows = num.pca, skip=1)

# plot screeplot
pdf(file = file.path(out.path, out.plotname))
  ggplot(pca, aes(seq_along(pca$V1), V1)) + 
    geom_line() + 
    ylab("Eigenvalue") +
    xlab("Principal component") +
    scale_x_continuous(breaks = seq_along(1:num.pca)) +
    ggtitle(paste(plot.header, "- screeplot - PCA final"))
dev.off()
