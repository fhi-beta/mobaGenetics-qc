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
input.sampleSheet = args[1]
output.rgset = args[2]

# ---- 1.Read in data
message("\nReading in sample sheet...\n")

arrays <- readRDS(input.sampleSheet)
# print(dim(arrays))

# Nice command to have in case:
# arrays$Basename = paste(plateDir, "/", paste(arrays$Slide, arrays$SentrixPosition, sep = "_"), sep = "")

message("Dimension of sample sheet: \n")
print(dim(arrays))

# ---- 2. Read data into RGChannelSets
message("\n")
message("Building RGChannelSet from array data...\n")

RGSet <- suppressWarnings(read.metharray.exp(targets=arrays, extended=FALSE))

# message("print size of session: \n")
# print(mem_used())

# ---- 4. Save the RGSet to disk
message("Saving the RGSet to disk...\n")

saveRDS(RGSet, file=output.rgset)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("Building RGSet Done!")

