# ---- 0. Load dependencies and parse arguments

# print start time of script:
start_time = Sys.time()
message(paste0("The script was started at: \n", start_time, "\n\n"))

# Parse arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

# activate renv if renv = TRUE
if(as.logical(args[1])){
	source("renv/activate.R")
}

args = args[-1]

message("Loading script dependencies...\n")

# Suppres package load messages
suppressMessages({
        library(minfi, quietly=TRUE)
	library(pryr, quietly=TRUE)
})


message("Printing arguments from snakemake...\n")
print(args)

# improve readability by unpacking args to variables:
input.rgset = args[1]
input.cross_hybridizing = args[2]

params.filter_cross_hybridizing = args[3]
params.filter_polymorphic_CpGs = args[4]
params.detection_p = args[5]

output.probes_removed = args[6]
output.SNPbetas = args[7]
output.rgset_filtered_probes = args[8]
output.pvalueMat = args[9]

# ---- 1.Read in RGSet
message("\nReading in RGSet...\n")

RGSet <- readRDS(input.rgset)

message("Dimension of RGset: \n")
print(dim(RGSet))

message("\n")
message(paste0("Filter out SNP probes: TRUE \nFilter out cross hybridization probes: ", params.filter_cross_hybridizing, "\nFilter out polymorphic CpGs: ", params.filter_polymorphic_CpGs, "\n"))

SNPbetas = getSnpBeta(RGSet)

cross_hybridizing = NULL
if(params.filter_cross_hybridizing){
	cross_hybridizing = read.table(input.cross_hybridizing)
	cross_hybridizing = as.character(cross_hybridizing[, 1])
}

polymorphic = NULL
if(params.filter_polymorphic_CpGs){
	snpDF = getSnpInfo(RGSet)
	wh = Reduce(union, lapply(c("CpG_maf", "SBE_maf"), function(xx){which(snpDF[, xx] >= 0.05)}))
	polymorphic = rownames(snpDF[sort(wh),])
}



RGSet = subsetByLoci(RGSet, excludeLoci = unique(c(cross_hybridizing, polymorphic)), keepControls = TRUE, keepSnps = FALSE)

pValMat = detectionP(RGSet)

message("\nDimension of RGset after filtering: \n")
print(dim(RGSet))

message("Dimension of pvalue-matrix: \n")
print(dim(pValMat))

message("\nCreate and store .csv of probes removed...")
probes_removed_mat = matrix("c", ncol = 2, nrow = length(unique(c(cross_hybridizing,polymorphic))))
rownames(probes_removed_mat) = unique(c(cross_hybridizing, polymorphic))
colnames(probes_removed_mat) = c("CpG", "Reason")
probes_removed_mat[,1] = rownames(probes_removed_mat)
probes_removed_mat[polymorphic, 2] = "polymorphic"
probes_removed_mat[cross_hybridizing, 2] = "cross_hybridizing"

write.csv(probes_removed_mat, file = output.probes_removed, row.names = F)

message("\nSave SNPbetas, filtered RGSet and pvalue-matrix...")

saveRDS(SNPbetas, file = output.SNPbetas)
saveRDS(RGSet, file = output.rgset_filtered_probes)
saveRDS(pValMat, file = output.pvalueMat)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("\nFiltering RGset probes done! \n\n")




