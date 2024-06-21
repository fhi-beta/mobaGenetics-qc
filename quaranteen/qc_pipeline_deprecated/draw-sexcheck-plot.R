#!/usr/bin/env Rscript
# This script generates plot images for SEXCHECK statistic.
# Args:
# 1. .sexcheck file from PLINK with Fstat
# 2. dir/ to store plots
# 3, 4. upper and lower thresholds for Fstat

library(ggplot2)

args = commandArgs(TRUE)
statfile = args[1]
plotdir = args[2]
thr_upp = as.numeric(args[3])
thr_low = as.numeric(args[4])

f_het_x = read.table(statfile, h=T)

p = ggplot(f_het_x, aes(x=F)) + geom_histogram() +
    stat_bin(aes(label=ifelse(..count..==0, "", ..count..), y=..count..),
             geom="text", vjust=-0.5, size=3) +
    geom_jitter(data=f_het_x[f_het_x$STATUS!="OK",], aes(y=1, x=F), color="red", width=0, height=10) +
    theme_bw() + ggtitle("Sample X chr. F-statistic") + xlab(NULL)
p = p + geom_vline(aes(xintercept=thr_upp), colour="red")
p = p + geom_vline(aes(xintercept=thr_low), colour="red")
ggsave(paste(plotdir, "fstat.png", sep="_"), width=16, height=10, units="cm", plot=p)

