args = commandArgs(trailingOnly=TRUE)

hetfile = args[1]
thr = args[2]
outfile = args[3]

het = read.table(hetfile, header = T)
het$HET_RATE=(het$"N.NM."-het$"O.HOM.")/het$"N.NM."
het_fail=subset(het,(het$HET_RATE>mean(het$HET_RATE)+as.numeric(thr)*sd(het$HET_RATE)))
het_fail$HET_DST=(het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)
write.table(het_fail[,1:2],file = paste0(outfile),row.names=F,quote=F)
