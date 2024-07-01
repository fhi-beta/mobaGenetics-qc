#!/usr/bin/Rscript

## converts long format table into wide table
# args: 1 - folder with exclusions_long file
library(dplyr)
library(tidyr)
options(stringsAsFactors = F)
args = commandArgs(TRUE)
indir = args[1]

# indir = "/mnt/HUNT/erc-genotypes-results/"
l = read.table(paste(indir, "exclusions_ind_long.txt", sep=""), h=T)
l$FID = as.character(l$FID)

# reorder stages, skip the temporary ones
## remove initial duplicate exclusion
l$STAGE = factor(l$STAGE, levels = c("superclean", "pca", "cleansex", "verify"))
l$FILTER = factor(l$FILTER, levels = c("MIND0.05", "INFERPED", "PCA", "SEXCHECK", 
                                       "MIND0.02", "HET4", "HET_RARE4",  "PIHAT_ACCUM", "PIHAT_THR"))
l = filter(l, !is.na(STAGE), !is.na(FILTER)) %>%
    arrange(STAGE, FILTER)

# read in all remaining individuals
f_all = read.table(paste(indir, "both/tmp/inferped_all.fam", sep=""))
f_all$V1 = as.character(f_all$V1)

# attach all-pass individuals, create PASS/FAIL column
print(head(f_all))
print(head(l))
f_pass = anti_join(f_all, l, by=c("V2"="IID"))[,c(1, 2)]
f_pass$STAGE = f_pass$FILTER = NA
f_pass$STATUS = "P"
l$STATUS = "F"
l = l[,-5]
colnames(f_pass) = colnames(l)
l = bind_rows(l, f_pass)

# make unique stage+filter column, convert to wide
w = unite(l, STAGE_FILTER, STAGE, FILTER, sep="__") %>%
    mutate(STAGE_FILTER=factor(STAGE_FILTER, levels=unique(STAGE_FILTER))) %>%
    unique() %>%
    spread(STAGE_FILTER, STATUS) %>%
    select(-matches("NA__NA"))

# if an excluded sample actually passed some tests, fill those with P
w = data.frame(t(apply(w, 1, function(x){
    m = match("F", x[1:11])
    if(!is.na(m)){
        x[3:(m)] = "P"
        x[m] = "F"
        if(m<11) x[(m+1):11] = NA
    } else x[3:11] = "P"
    x
})))

# add flag to mark reliable genotypes
w$genotypesOK = w$superclean__MIND0.05=="P"
w$coreLMM = !apply(w, 1, function(x) any(x[3:9]=="F", na.rm=T))
w$coreUNRELATED = !apply(w, 1, function(x) any(x[3:11]=="F", na.rm=T))

# read in pedigree flags
p = read.table(paste(indir, "both/rep/phenotypesOK.txt", sep=""), h=T)

wp = left_join(w, p, by="IID")

## TODO: use input recode files

# read in files to identify samples as founder/offspring
id_offspring = read.table("/media/local-disk2/helgeland/rotterdam1/recode-files/recode-parents-rotterdam1.txt", h=F)

wp$ROLE = ifelse(wp$IID %in% id_offspring[,2], "OFFSPRING", "FOUNDER")

write.table(wp, paste(indir, "sample_flag_list.txt", sep=""),
            sep="\t", col.names=T, row.names=F, quote=F)
