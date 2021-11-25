# ---- 0. Load dependencies and parse arguments

message("Loading script dependencies...\n")

# Suppres package load messages
suppressMessages({
        library(minfi, quietly=TRUE)
	library(pryr, quietly=TRUE)
})

# Parse arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 
# message("Print arguments from snakemake: \n")
# print(args)

message("Printing arguments from snakemake...\n")
print(args)

# improve readability by unpacking args to variables:
input.rgset = args[1]
input.cross_hybridizing = args[2]

params.filter_cross_hybridizing = args[3]
params.filter_polymorphic_CpGs = args[4]
params.detection_p = args[5]
params.frac_poor_quality = args[6]

output.probes_removed = args[7]
output.SNPbetas = args[8]
output.rgset_filtered_probes = args[9]

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

message(paste0("\n", "Filter out CpGs with detection p-value less than ", params.detection_p, " across more than ", params.frac_poor_quality, " of the samples \n"))

pValMat = detectionP(RGSet)
drop = rowSums(pValMat > as.numeric(params.detection_p)) > ncol(RGSet)*as.numeric(params.frac_poor_quality)
message("Amount of droped probes based on detection p val: ")
print(table(drop))
filter_out = names(drop[drop == TRUE])

RGSet = subsetByLoci(RGSet, excludeLoci = unique(c(cross_hybridizing, polymorphic, filter_out)), keepControls = TRUE, keepSnps = FALSE)

message("\nDimension of RGset after filtering: \n")
print(dim(RGSet))

message("\nCreate and store .csv of probes removed...")
probes_removed_mat = matrix("c", ncol = 2, nrow = length(unique(c(cross_hybridizing,polymorphic,filter_out))))
rownames(probes_removed_mat) = unique(c(cross_hybridizing, polymorphic, filter_out))
colnames(probes_removed_mat) = c("CpG", "Reason")
probes_removed_mat[,1] = rownames(probes_removed_mat)
probes_removed_mat[filter_out, 2] = "high_detection_p"
probes_removed_mat[polymorphic, 2] = "polymorphic"
probes_removed_mat[cross_hybridizing, 2] = "cross_hybridizing"

write.csv(probes_removed_mat, file = output.probes_removed, row.names = F)

message("\nSave SNPbetas and filtered RGSet...")

saveRDS(SNPbetas, file = output.SNPbetas)
saveRDS(RGSet, file = output.rgset_filtered_probes)

message("\nFiltering RGset probes done! \n\n")




