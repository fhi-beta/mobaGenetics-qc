# ---- 0. Parse Snakemake arguments

# print start time of script:
start_time = Sys.time()
message(paste0("The script was started at: \n", start_time, "\n\n"))

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
	library(ggplot2, quietly=TRUE)
})


message("Printing arguments from snakemake...\n")
print(args)
message('\n')

# unpack args to improve readability:
input.rgset = args[1]
input.methset = args[2]
input.bmiqed = args[3]
output.plot = args[4]

# load the different data sets we need for plotting Raw values, Noob values and BMIQed values
bmiqed = readRDS(input.bmiqed)
set.seed(10)
cpg_sub = sample(rownames(bmiqed), size = 50000)
bmiqed = bmiqed[cpg_sub,]
invisible(gc())
message("\nDimension of BMIQed_data:")
print(dim(bmiqed))

rgset = readRDS(input.rgset)
rgset = getBeta(rgset)
rgset = rgset[cpg_sub,]
invisible(gc())
message("\nDimension of RGset:")
print(dim(rgset))

methset = readRDS(input.methset)
methset = getBeta(methset)
methset = methset[cpg_sub,]
invisible(gc())
message("\nDimension of Methylset:")
print(dim(methset))


overlap_samp = intersect(intersect(colnames(rgset), colnames(methset)), colnames(bmiqed))
# in case some samples got filtered out in manual QC step

# calculate the differences for each position between the samples, and take the mean of it:
Noob_Raw_mean = rowMeans(rgset[,overlap_samp] - methset[,overlap_samp], na.rm = TRUE)
Noob_Raw_sd = rowSds(rgset[,overlap_samp] - methset[,overlap_samp], na.rm = TRUE)

Bmiq_Noob_mean = rowMeans(methset[,overlap_samp] - bmiqed[,overlap_samp], na.rm = TRUE)
Bmiq_Noob_sd = rowSds(methset[,overlap_samp] - bmiqed[,overlap_samp], na.rm = TRUE)


message('Saving plots to file...\n')

pdf(output.plot, height=6, width=8.5)
ggplot(data.frame(x = Noob_Raw_mean), aes(x = x)) + geom_histogram(bins = 50) + ggtitle("Noob beta values vs raw") + xlab("Mean of methylation differences between samples")
ggplot(data.frame(x = Noob_Raw_sd), aes(x = x)) + geom_histogram(bins = 50) + ggtitle("Noob beta values vs raw") + xlab("std. dev. of methylation differences between samples")
ggplot(data.frame(x = Bmiq_Noob_mean), aes(x = x)) + geom_histogram(bins = 50) + ggtitle("Bmiq beta values vs Noob") + xlab("Mean of methylation differences between samples")
ggplot(data.frame(x = Bmiq_Noob_sd), aes(x = x)) + geom_histogram(bins = 50) + ggtitle("Bmiq beta values vs Noob") + xlab("std. dev. of methylation differences between samples")
dev.off()

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time


message("Plotting done!\n\n")


