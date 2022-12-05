#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

pcainput = args[1]
outfile = args[2]

# Read in the PCA input from Eigenstrat. The format expects FID and IID to be
# separated by ":". 
pca <- read.table(pcainput, skip = 1, header = F)
pca <- separate(pca, V1, into = c("V0","V1"), sep = ":", extra = "merge")

# Plot
p <- ggplot(data= pca, aes(V2,V3)) + 
    geom_point() + 
    labs(x="PC1", y="PC2")
ggsave(outfile, width=10, height=10, units="cm", plot=p)
