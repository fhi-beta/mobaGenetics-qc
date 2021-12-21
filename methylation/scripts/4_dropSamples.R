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
    library(IlluminaHumanMethylationEPICmanifest, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
})

message("Printing arguments from snakemake...\n")
print(args)
# improve readability by unpacking args to variables:
input.rgset = args[1]
input.NA_control_probes = args[2]

params.detection_pvalue = args[3]
params.bisulphite_conversion_rate = args[4]
params.poor_quality = 1 - as.numeric(args[5])

output.detectionPvalues = args[6]
output.probeQC = args[7]
output.sampleQC = args[8]
output.bisulphiteQC = args[9]
output.rgset_filtered = args[10]
output.removed_samples = args[11]

# ---- 1. Load RGChannelSet
message("\nReading in RGChannelSet from: ", input.rgset, '...\n')
rgSet <- readRDS(input.rgset)


# ---- 2. Calculate detection p-value, proportion/count failed probes per sample, probes with proportion failed samples 
message("Calculating QC metrics...\n")
source('scripts/functions/summarizeDetectionPvalueQC.R')

# make initial filtering of samples data frame for storing of samples that are removed, and why
# dropped_samples_df <- data.frame(matrix(ncol=2,nrow=0))
# colnames(dropped_samples_df) <- c('sample','reason')
dropped_samples_df = read.csv(file = input.NA_control_probes)
if(nrow(dropped_samples_df) > 0){
	rgSet = rgSet[,setdiff(colnames(rgSet), dropped_samples_df[,1])]
	print(dim(rgSet))
	message(paste0("Filtered out ", nrow(dropped_samples_df), " samples due to NA on control probes... \n"))
}

# Remove failed samples (too low concentration of DNA) and illumina controls, as these are not used
# only relevant for met001 and met002, which has a "Comment" column with this information
pheno = pData(rgSet)
if(sum(colnames(pheno) == "Comment")){
	message(paste0("Dimension of rgSet before intial filtering: ", dim(rgSet)[1], " ", dim(rgSet)[2]))
	comments = tolower(pheno[, "Comment"])
	comments = gsub(" ", "", comments) # remove spaces
	names(comments) = rownames(pheno)
	exclude = names(comments[grep("illuminactrl", comments)])
	reason = rep("illumina_control", times = length(exclude))
	exclude = c(exclude, names(comments[grep("failedsample", comments)]))
	reason = c(reason, rep("failed_sample", times = length(comments[grep("failedsample", comments)])))

	dropped_samples_df = rbind(dropped_samples_df, data.frame(sample=exclude, reason=reason))
	message(paste0("Remove ", length(exclude), " failed samples(too low concentration) and illumina controls.."))
	rgSet = rgSet[, setdiff(colnames(rgSet), exclude)]
	message(paste0("Dimension of rgSet after intial filtering: ", dim(rgSet)[1], " ", dim(rgSet)[2]))
}

detect_p_val = as.numeric(params.detection_pvalue)
detectionQCResults <- summarizeDetectionPvalueQC(rgSet, pValue=detect_p_val)


# ---- 4. Construct result names and save QC results to disk

detectionPOut <- output.detectionPvalues
sampleQCOut <- output.sampleQC
probeQCOut <- output.probeQC

paths <- c(detectionPOut, sampleQCOut, probeQCOut)

for (i in seq_along(detectionQCResults)) {
    message(paste0('\tSaving QC results to: ', paths[i], '...\n'))
    fwrite(detectionQCResults[[i]], file=paths[i])
}


# ---- 5. Filter failed samples from RGChannelSet and save to disk
message(paste0('Filtering RGChannelSet to samples with >90% CpGs detected at alpha ', detect_p_val, '...\n'))
# Drop rownames column so all values are numeric and we can use colMeans
detPvals <- detectionQCResults$detectionPvals[, -'rownames']
detected <- detPvals < detect_p_val


numSamplesPassedQC <- sum(colMeans(detected) > params.poor_quality)
message(paste0(numSamplesPassedQC, ' out of ', ncol(detPvals), ' samples passed detection p-value filtering  at alpha of ', detect_p_val, '...\n'))

failedSamples <- which(colMeans(detected) < params.poor_quality)

keepSamples <- colMeans(detected) > params.poor_quality

rgSet <- rgSet[, keepSamples]
invisible(gc())

df <- data.frame(sample=names(failedSamples), reason=rep('pvalue', length=length(failedSamples)))
dropped_samples_df <- rbind(dropped_samples_df, df)


# ---- 6. Calculate the bisulphite conversion rates per sample and save to disk
message("Calculating sample bisulphite conversion rates...\n")
bisulphiteConversion <- bscon(rgSet)
bisulphiteConversionDT <- data.table(
                            'sample'=names(bisulphiteConversion),
                            'conversion_rate'=bisulphiteConversion
                            )


fwrite(bisulphiteConversionDT, file=output.bisulphiteQC)

# ---- 7. Determine which how many and which samples fail bisulphite conversion qc
conversion_minimum <- as.numeric(params.bisulphite_conversion_rate)

numPassed <- sum(bisulphiteConversionDT$conversion_rate > conversion_minimum)

message(paste0(numPassed, ' out of ', nrow(bisulphiteConversionDT),
               ' samples passed qc at bisulphite conversion rate of ',
               conversion_minimum, '%...\n'))

# ---- 8. Drop samples failing the bisulphite conversion cut-off
message("Subsetting RGChannelSet to samples passing bisulphite QC...\n")
keepSamples <- bisulphiteConversion > conversion_minimum

rgSet <- rgSet[, keepSamples]
invisible(gc())
message('\n')


# ---- 9. Save RGChannelSet to disk
message(paste0("Saving RGChannel set to ", output.rgset_filtered , '...\n'))
saveRDS(rgSet, file=output.rgset_filtered)

dropped_samples_df

failedSamples2 <- which(bisulphiteConversion <= conversion_minimum)
df <- data.frame(sample=names(failedSamples2), reason=rep('bisulphite', length=length(failedSamples2)))
dropped_samples_df <- rbind(dropped_samples_df, df)

# ---- 10. Save dropped samples to disk
message(paste0("Saving dropped samples to ", output.removed_samples, '...\n'))
write.csv(dropped_samples_df, file = output.removed_samples, row.names = F)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("Filtering RGSet samples done!\n\n")
