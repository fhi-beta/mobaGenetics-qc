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
})

message("Printing arguments from snakemake...")
print(args)
message("\n")

# unpack args to improve readability:
input.ratioset = args[1]
output.filteredRatioset = args[2]
output.addedPheno = args[3]

# ---- 1. Read in data
message(paste0('\nReading in GRset: ', input.ratioset, '...\n'))
GRset = readRDS(input.ratioset)

sex_prediction = getSex(GRset)

pheno = pData(GRset)
if(sum(colnames(pheno) == "Sex")){
	if("FALSE" %in% as.character(pheno$Sex)){
	pheno$Sex = substr(as.character(pheno$Sex), 1, 1)
	}
}

pheno_added = cbind(pheno, sex_prediction[rownames(pheno), ])
message("\nDimension of GRset before filtering: \n")
print(dim(GRset))


# if the sample sheet has a column with specified "Sex", check for discordancy
if(sum(colnames(pheno) == "Sex")){
	message(paste0("There are ", sum(is.na(pheno$Sex)), " NAs in the Sex column.."))
	if(sum(is.na(pheno$Sex)) > 0){
	}else{
	discordant = sum(pheno_added[, "Sex"] != pheno_added[, "predictedSex"])
	message(paste0("\nThere are ", discordant, " miss-matches between Sex and predicted sex\n"))
	message("Filter out sex discordant samples...")
	keep = rownames(pheno_added[pheno_added[, "Sex"] == pheno_added[, "predictedSex"],])
	GRset = GRset[, keep]
	}
}
message("Dimension of GRset after filtering: \n")
print(dim(GRset))

message("\nSave filtered GRset and pheno data with added sex prediction column...")

saveRDS(GRset, file = output.filteredRatioset)
write.csv(pheno_added, file = output.addedPheno, row.names = T)

message("\nPredicting sex and checking for discordancy done!\n")
