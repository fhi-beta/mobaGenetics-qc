#!/usr/bin/env Rscript
# This script generates plot imagess from PLINK's summary stat files.
# Args:
# 1. summary stat filestem
# 2. stage name
# 3-7. Thresholds for LMISS, IMISS, HET, log(HWE) and MAF filters

library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = TRUE)
statfiles = args[1]
plotdir = args[2]
stage = args[3]
thr_lmiss = as.numeric(args[4])
thr_imiss = as.numeric(args[5])
thr_het = as.numeric(args[6])
thr_hwe = -as.numeric(args[7])
thr_maf = as.numeric(args[8])

zoomtheme = theme(axis.title=element_blank(), title=element_blank(),
                  plot.margin = grid::unit(c(3, 3, 3, 3), "mm"),
                  panel.background = element_rect(color="grey50", fill="white"))

# snp missingness plot
if(!is.na(thr_lmiss)){
	f_lmiss = read.table(paste(statfiles, ".lmiss", sep=""), header = T)
	f_lmiss = f_lmiss[order(1-f_lmiss$F_MISS),]
	f_lmiss$markers = seq_along(f_lmiss$SNP)
	p = ggplot(f_lmiss, aes(x=markers, y=1-F_MISS)) + geom_line() + geom_point(shape=18, size=1) +
		ggtitle("Marker missingness") + 
		geom_hline(aes(yintercept=1-thr_lmiss), colour="red")
	pzoom = p + coord_cartesian(xlim = c(0, 0.05*nrow(f_lmiss)), ylim = c(0.9, 1)) + zoomtheme
	limlow = min(1-thr_lmiss, 1-max(f_lmiss$F_MISS, na.rm=T))
	pall = p + annotation_custom(ggplotGrob(pzoom),
	                             xmin=0.5*nrow(f_lmiss), xmax=nrow(f_lmiss),
			  ymin=limlow, ymax=0.5*(limlow+1))
	ggsave(paste(plotdir, stage, "lmiss.png", sep="_"), width=16, height=10, units="cm", plot=pall)
}

# sample missingness plot
if(!is.na(thr_imiss)){
	f_imiss = read.table(paste(statfiles, ".imiss", sep=""), header = T)
	f_imiss = f_imiss[order(1-f_imiss$F_MISS),]
	f_imiss$samples = seq_along(f_imiss$IID)
	p = ggplot(f_imiss, aes(x=samples, y=1-F_MISS)) + geom_line() + geom_point(shape=18, size=1) +
	    ggtitle("Sample missingness") +
	    geom_hline(aes(yintercept=1-thr_imiss), colour="red")
	pzoom = p + coord_cartesian(xlim = c(0, 0.05*nrow(f_imiss)), ylim = c(0.9, 1)) + zoomtheme
	limlow = min(1-thr_imiss, 1-max(f_imiss$F_MISS, na.rm=T))
	pall = p + annotation_custom(ggplotGrob(pzoom),
	                             xmin=0.5*nrow(f_imiss), xmax=nrow(f_imiss),
			  ymin=limlow, ymax=0.5*(limlow+1))
	ggsave(paste(plotdir, stage, "imiss.png", sep="_"), width=16, height=10, units="cm", plot=pall)
}

# heterozygosity plots
if(!is.na(thr_het)){
	f_het = read.table(paste(statfiles, ".het", sep=""), header = T)
	f_het_rare = read.table(paste(statfiles, "_rare.het", sep=""), header = T)
	f_het$RATE = (f_het$N.NM.-f_het$O.HOM.)/f_het$N.NM.
	thr_het_u = mean(f_het$RATE) + thr_het*sd(f_het$RATE)
	p = ggplot(f_het, aes(x=RATE)) + geom_histogram() +
	    stat_bin(aes(label=ifelse(..count..==0, "", ..count..), y=..count..),
	             geom="text", vjust=-0.5, size=3) +
	    theme_bw() + ggtitle("Sample heterozygosity") + xlab(NULL) +
	    geom_vline(aes(xintercept=thr_het_u), colour="red")
	ggsave(paste(plotdir, stage, "het.png", sep="_"), width=16, height=10, units="cm", plot=p)
	
	f_het_rare$RATE = (f_het_rare$N.NM.-f_het_rare$O.HOM.)/f_het_rare$N.NM.
	thr_het_rare_u = mean(f_het_rare$RATE) + thr_het*sd(f_het_rare$RATE)
	p = ggplot(f_het_rare, aes(x=RATE)) + geom_histogram() +
	    stat_bin(aes(label=ifelse(..count..==0, "", ..count..), y=..count..),
	             geom="text", vjust=-0.5, size=3) +
	    theme_bw() + ggtitle("Sample heterozygosity") + xlab(NULL) + 
	    geom_vline(aes(xintercept=thr_het_rare_u), colour="red")
	ggsave(paste(plotdir, stage, "het_rare.png", sep="_"), width=16, height=10, units="cm", plot=p)
}

if(!is.na(thr_hwe)){
	f_hwe = read.table(paste(statfiles, ".hwe", sep=""), header = T)
	f_hwe = f_hwe[which(f_hwe$P!=0),]
	f_hwe$logP_observed = sort(-log(f_hwe$P, 10))
	f_hwe$logP_expected = sort(-log(runif(nrow(f_hwe)), 10))
	p = ggplot(f_hwe) + geom_point(aes(x=logP_expected, y=logP_observed)) +
	    geom_abline(slope=1, intercept=0, col="black") +
	    ggtitle("Q-Q plot of Hardy-Weinberg equilibrium P-values") + 
	    geom_hline(aes(yintercept=thr_hwe), colour="red")
	pzoom = ggplot(f_hwe, aes(x=logP_observed)) + geom_histogram() + coord_flip() +
	    theme_minimal() + xlab(NULL) + ggtitle("") +
	    stat_bin(aes(label=ifelse(..count..==0 | ..count..>1000, "", ..count..), y=..count..),
	                geom="text", hjust=-0.3, size=3) +
	    scale_x_continuous(breaks = scales::pretty_breaks(n=3))
	pall = grid.arrange(p, pzoom, ncol=2, widths=c(0.7, 0.3))
	ggsave(paste(plotdir, stage, "hwe.png", sep="_"), width=16, height=10, units="cm", plot=pall)
}

# f_maf = read.table(paste(statfiles, ".frq", sep=""), header = T)
# p = ggplot(f_maf) + geom_histogram(aes(x=MAF), binwidth=0.02) +
#     theme_bw() + ggtitle("SNP minor allele frequency") + xlab(NULL)
# ggsave(paste(plotdir, stage, "maf.png", sep="_"), width=16, height=10, units="cm", plot=p)

