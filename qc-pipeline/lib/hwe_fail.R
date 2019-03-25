#!/usr/bin/Rscript

# Get arguments from command line
args = commandArgs(trailingOnly=TRUE)

# PLINK HWE-statistics 
hardyfile = args[1]

# Threshold for HWE P-value (markers with P-value below this threshold will be excluded)
thr = as.numeric(args[2])

# Path to output file containing only markers below P-value
outfile = args[3]

# If not NA - path to output plot
plotfile = args[4]

# Read PLINK HWE-table
hardy <- read.table(hardyfile, header = T, stringsAsFactors=F)

# Subset only markers below P-value
s <- subset(hardy, P < thr)

# If no markers are in the subset, only output an empty file
# else write only the markers to output
if(nrow(s)==0) {  file.create(outfile)} else write.table(s$SNP,file = outfile, row.names=F, col.names = F, quote=F)

if(plotfile!="NA"){
	library(ggplot2)
	library(gridExtra)

	hardy = hardy[which(hardy$P!=0),]
	hardy$logP_observed = sort(-log(hardy$P, 10))
	hardy$logP_expected = sort(-log(runif(nrow(hardy)), 10))
	p = ggplot(hardy) + geom_point(aes(x=logP_expected, y=logP_observed)) +
		geom_abline(slope=1, intercept=0, col="black") +
		ggtitle("Q-Q plot of Hardy-Weinberg equilibrium P-values, X chr.") +
		geom_hline(aes(yintercept = -log(thr, 10)), colour="red")
	pzoom = ggplot(hardy, aes(x=logP_observed)) + geom_histogram() + coord_flip() +
	    theme_minimal() + xlab(NULL) +
	    stat_bin(aes(label=ifelse(..count..==0 | ..count..>200, "", ..count..), y=..count..),
                geom="text", hjust=-0.3, size=3)
	pall = grid.arrange(p, pzoom, ncol=2, widths=c(0.7, 0.3))
	ggsave(plotfile, width=16, height=10, units="cm", plot=pall)
}

