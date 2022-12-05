#!/usr/bin/env Rscript

# Parse arguments
args <- commandArgs(TRUE)

snptable.path                     <- args[1]

threshold.cluster.sep             <- args[2]
cluster.sep.exclude.outlist.path  <- args[3]

threshold.tenperc.gc.score        <- args[4]
tenperc.gc.exclude.outlist.path   <- args[5]

threshold.aa.t.dev                <- args[6]
aa.t.dev.exclude.outlist.path     <- args[7]

# Load libraries
library(data.table)

# Load data 

# load snp table
snptable <- fread(snptable.path, header = T, stringsAsFactors = F)

# make names R safe for convenience
names(snptable) <- make.names(names(snptable))

# CLUSTER SEPARATION
cluster.sep.col <- names(snptable)[grep(pattern = 'bpm.Cluster.Sep', names(snptable))]
cl.sep.exclude <- subset(snptable, get(cluster.sep.col) < threshold.cluster.sep)
write.table(cl.sep.exclude$Name, file = cluster.sep.exclude.outlist.path, col.names = F, row.names = F, quote = F)

# 10 % GC SCORE
tenperc.gc.col <- names(snptable)[grep(pattern = 'X10..GC', names(snptable))]
tenperc.gc.exclude <- subset(snptable, get(tenperc.gc.col) < threshold.tenperc.gc.score)
write.table(tenperc.gc.exclude$Name, file = tenperc.gc.exclude.outlist.path, col.names = F, row.names = F, quote = F) 

# AA THETA DEV
aa.theta.col <- names(snptable)[grep(pattern = 'bpm.AA.T.Dev', names(snptable))]
aa.theta.dev.exclude <- subset(snptable, get(aa.theta.col) > threshold.aa.t.dev)
write.table(aa.theta.dev.exclude$Name, file = aa.t.dev.exclude.outlist.path, col.names = F, row.names = F, quote = F)

nrow(subset(snptable, get(aa.theta.col) > 0.025))

#ggplot(snptable) + geom_histogram(aes(HumanCoreExome.24v1.0_A.bpm.Cluster.Sep), bins=100) + ggtitle("Cluster separation all markers")
#ggplot(snptable) + geom_histogram(aes(HumanCoreExome.24v1.0_A.bpm.Cluster.Sep), bins=100) + xlim(c(0, 0.99)) + 
#  ggtitle("Cluster separation for markers with cluster sep 0 - 0.99") + geom_vline(xintercept= 0.4, col='red')

#ggplot(snptable) + geom_histogram(aes(X10..GC), bins=200)

#ggplot(subset(snptable, Call.Freq > 0.98)) + geom_point(aes(Call.Freq, X10..GC))

#ggplot(arrange(snptable, Call.Freq)) + geom_point(aes(seq_along(Call.Freq), Call.Freq)) + ggtitle("Call frequency") + geom_hline(yintercept= 0.98, col='red')

