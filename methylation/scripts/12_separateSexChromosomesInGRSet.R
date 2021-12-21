# ---- 0. Load dependencies

# print start time of script:
start_time = Sys.time()
message(paste0("The script was started at: \n", start_time, "\n\n"))

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
input.pvalueMat = args[2]
output.ratiosetNoSex = args[3]
output.pvalMatNoSex = args[4]
output.pvalMatNoSexRds = args[5]
output.betaMatrixOnlySex = args[6]
output.betaMatrixSexChrRds = args[7]
output.pvalMatOnlySex = args[8]
output.pvalMatOnlySexRds = args[9]

# ---- 1. Read in data
message(paste0('Reading in GRSet from: ', input.ratioset, '...\n'))
grSet <- readRDS(input.ratioset)
pvalMat = readRDS(input.pvalueMat)

# ---- 2. Separate sex chromosomes and save to GenomicRatioSet
message(paste0('Separating sex chromosomes from GenomicRatioSet and saving to: ', output.ratiosetNoSex, ' and ', output.betaMatrixOnlySex, '\n'))
grSet_annotation <- getAnnotation(grSet)
sex_chr <- (grSet_annotation$chr %in% c('chrX', 'chrY'))

message('Sex chromosome: ')
print(table(sex_chr))

grSet_no_sex_chr <- grSet[!sex_chr, ]

message("Dimension of autosomal data frame: \n")
print(dim(grSet_no_sex_chr))
message("Dimension of autosomal pvalue matrix: \n")
print(dim(pvalMat[rownames(grSet_no_sex_chr), colnames(grSet_no_sex_chr)]))

message("Save autosomal CpGs for further QC analysis...\n")
saveRDS(grSet_no_sex_chr, output.ratiosetNoSex)

message("Save pvalue matrix for autosomal CpGs and sex chromsome CpGs separately...\n")
write.csv(pvalMat[rownames(grSet_no_sex_chr), colnames(grSet_no_sex_chr)], file = output.pvalMatNoSex, row.names = T)
saveRDS(pvalMat[rownames(grSet_no_sex_chr), colnames(grSet_no_sex_chr)], file = output.pvalMatNoSexRds)


# store the beta matrix directly, as no more qc will be done on the sex chromosomes
beta_matrix_sex_chr <- getBeta(grSet[sex_chr, ])
write.csv(beta_matrix_sex_chr, file = output.betaMatrixOnlySex, row.names = T)
saveRDS(beta_matrix_sex_chr, file = output.betaMatrixSexChrRds)

message("Dimension of sex chromosome specific CpG data frame: \n")
print(dim(beta_matrix_sex_chr))
message("Dimension of sex chromosome specific pvalue matrix: \n")
print(dim(pvalMat[rownames(beta_matrix_sex_chr), colnames(beta_matrix_sex_chr)]))

write.csv(pvalMat[rownames(beta_matrix_sex_chr), colnames(beta_matrix_sex_chr)], file = output.pvalMatOnlySex, row.names = T)
saveRDS(pvalMat[rownames(beta_matrix_sex_chr), colnames(beta_matrix_sex_chr)], file = output.pvalMatOnlySexRds)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("\nSeparating sex chromosomes done!\n\n")

