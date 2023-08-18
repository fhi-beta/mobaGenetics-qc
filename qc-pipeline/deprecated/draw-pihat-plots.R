#!/usr/bin/env Rscript
# This script makes a Z1-Z0 summary plot from PIHAT statistics

# Args:
# 1. input .genome file
# 2. output .png file

# options(stringsAsFactors=F)
args = commandArgs(TRUE)
genofile = args[1]
outfile = args[2]
colorize = args[3]

library(ggplot2)

gf = read.table(genofile, h=T)
gf = unique(gf)
if(colorize=="TRUE"){
	p = ggplot(gf, aes(x=Z0, y=Z1, color=RT)) + geom_point(size=0.5) +
	    theme_bw() + guides(color=guide_legend(override.aes=list(size=5)))
} else {
	p = ggplot(gf, aes(x=Z0, y=Z1)) + geom_point(size=0.5) + theme_bw()
}
ggsave(outfile, width=16, height=10, units="cm", plot=p)

