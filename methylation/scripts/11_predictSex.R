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
})

message("Printing arguments from snakemake...")
print(args)
message("\n")

# unpack args to improve readability:
input.ratioset = args[1]
output.filteredRatioset = args[2]
output.addedPheno = args[3]
output.sex_plot = args[4]

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

message("Store plot of xMed and yMed with color on sex.. \n")
pdf(output.sex_plot, height=6, width=8.5)
plot(x=pheno_added$xMed, y=pheno_added$yMed, col = factor(pheno_added$Sex), ylab = "Y chromosome median", xlab = "X chromosome median")
dev.off()

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

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("\nPredicting sex and checking for discordancy done!\n")
