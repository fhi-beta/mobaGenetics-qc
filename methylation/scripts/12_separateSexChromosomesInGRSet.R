# ---- 0. Load dependencies

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality

# activate renv if renv = TRUE
if(as.logical(args[1])){
        source("renv/activate.R")
}

args = args[-1]

message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
})

message("Printing arguments from snakemake...\n")
print(args)
message("\n")

# unpack args to improve readability:
input.ratioset = args[1]
output.ratiosetNoSex = args[2]
output.betaMatrixOnlySex = args[3]
output.betaMatrixSexChrRds = args[4]

# ---- 1. Read in data
message(paste0('Reading in GRSet from: ', input.ratioset, '...\n'))
grSet <- readRDS(input.ratioset)

# ---- 2. Separate sex chromosomes and save to GenomicRatioSet
message(paste0('Separating sex chromosomes from GenomicRatioSet and saving to: ', output.ratiosetNoSex, ' and ', output.betaMatrixOnlySex, '\n'))
grSet_annotation <- getAnnotation(grSet)
sex_chr <- (grSet_annotation$chr %in% c('chrX', 'chrY'))

message('Sex chromosome: ')
print(table(sex_chr))

grSet_no_sex_chr <- grSet[!sex_chr, ]
saveRDS(grSet_no_sex_chr, output.ratiosetNoSex)

# store the beta matrix directly, as no more qc will be done on the sex chromosomes
beta_matrix_sex_chr <- getBeta(grSet[sex_chr, ])
write.csv(beta_matrix_sex_chr, file = output.betaMatrixOnlySex, row.names = T)
saveRDS(beta_matrix_sex_chr, file = output.betaMatrixSexChrRds)

message("\nSeparating sex chromosomes done!\n\n")

