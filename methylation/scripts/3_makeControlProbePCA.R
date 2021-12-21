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
	library(irlba, quietly=TRUE)
	library(ggplot2, quietly=TRUE)
})


message("Printing arguments from snakemake...\n")
print(args)


# improve readability by unpacking args to variables:

input.rgset = args[1]
output.PCAPlot = args[2]
output.controlProbePCA = args[3]
output.NA_control_probes = args[4]
output.boxplot_control_probes = args[5]
# load RGset
message("\nReading in RGSet...\n")

RGset <- readRDS(input.rgset)

message("Subset on control probes...\n")
# get probe info annotation for the control probes
probeInfo_control = getProbeInfo(RGset, type = "Control")

# get the data; intensity measures from red and green channel
r_intensity = getRed(RGset)
g_intensity = getGreen(RGset)

# subset on control probes
r_control = r_intensity[probeInfo_control$Address,]
g_control = g_intensity[probeInfo_control$Address,]

# log values to normalize, prior to PCA
r_control_log2 = log2(r_control)
g_control_log2 = log2(g_control)
comb_control = rbind(r_control_log2, g_control_log2)
comb_control[is.infinite(comb_control)] = NA  # make infinite values equal to NA, for easy removement

amount_NA_control_probes = colSums(is.na(comb_control)) > 0.25*nrow(comb_control)
filter_mat = matrix("c", nrow = sum(amount_NA_control_probes), ncol = 2)
colnames(filter_mat) = c("sample", "reason")
if(sum(amount_NA_control_probes) > 0){
	filter_samples = names(amount_NA_control_probes[amount_NA_control_probes])
	message(paste0("\n", length(filter_samples), 
		       " of the samples have more than 25% missing values for control probes. 
		       They will be filtered out in next rule \n"))

	filter_mat[,1] = filter_samples
	filter_mat[,2] = "NA_control_probes"
	comb_control = comb_control[,-which(colnames(comb_control) %in% filter_samples)]
}

# save the sample names and reason for filtering for the filtered samples
write.csv(filter_mat, file = output.NA_control_probes, row.names = F)

# remove NA values and infinite values 
comb_control = na.omit(comb_control)
message(paste0("Dimension of na.omit of combined red and green control probes: ", 
	       dim(comb_control)[1], ", ", dim(comb_control)[2]), "\n")

message("Do PCA...\n")
PCA = irlba(t(comb_control), nv = 2, center = colMeans(t(comb_control)), 
	    scale = colSds(t(comb_control)))

rownames(PCA$u) = colnames(comb_control)

message("Create plot and store PCA...\n")

pheno = pData(RGset)
colvector = pheno$Slide
names(colvector) = rownames(pheno)

pdf(output.PCAPlot)
plot(PCA$u, col = factor(colvector[rownames(PCA$u)]), ylab = "2. PC", xlab = "1. PC")
dev.off()

saveRDS(PCA, file = output.controlProbePCA)

message("\n")
message("Create boxplot of control probes..\n")
pdf(output.boxplot_control_probes, height = dim(comb_control)[2]/2)
ggplot(data = data.frame(methylation = as.vector(comb_control), 
			 sample = factor(rep(colnames(comb_control), 
					     each = dim(comb_control)[1]))), 
       aes(x = methylation, y = sample, group = sample)) + geom_boxplot()
dev.off()

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time


message("\nMaking control probes plots  one!!!\n\n")
