# ---- 0. Load dependencies and parse arguments
message("Load .Rprofile for set up... \n")
source(".Rprofile")

message("Loading script dependencies...\n")

# Suppres package load messages
suppressMessages({
        library(minfi, quietly=TRUE)
})

# Parse arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

message("printing arguments from snakemake...\n")
print(args)
message("\n")
# Unpack args to improve readability:
input.methset = args[1]
input.filtered_samples = args[2]
output.qcPlot = args[3]
output.filtered_methset = args[4]
output.filtered_sample_list = args[5]
# ---- 1.Read in MethylSet
message("Reading in MSet...\n")
methylSet <- readRDS(input.methset)


# ---- 2.Make and save QC plot
message("Making and saving QC plot...\n")
qc <- getQC(methylSet)
df = read.csv(file = input.filtered_samples)

filter = rownames(qc[(qc[,1] + qc[,2])/2 < 10.5,])
if(length(filter) > 0){
	df = rbind(df, data.frame(sample=filter, reason=rep("low_median_meth_or_unmeth_channel", times = length(filter))))
}

write.csv(df, file=output.filtered_sample_list, row.names=F)

methylSet = methylSet[,setdiff(colnames(methylSet), filter)]

saveRDS(methylSet, file = output.filtered_methset)


pdf(output.qcPlot, height=11, width=8.5)
plotQC(qc)
dev.off()

message("Making QC plot done!")	
