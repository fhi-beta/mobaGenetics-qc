# ---- 0. Load dependencies and parse arguments

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
input.plateDir = args[1]
output.rgset = args[2]

# ---- 1.Read in data
message("\nReading in plate data...\n")

# Notes: if unable to load all idats due to memory limitation

plateDir <- input.plateDir
arrays <- read.metharray.sheet(plateDir)
# Notes: if unable to load all idats due to memory limitation
# subset the arrays before loading the idats

message("Dimension of sample sheet: \n")
print(dim(arrays))

# test if the sample composition is of different sampleType (tissues)
if(sum(colnames(arrays) == "SampleType")){
        uniq = unique(arrays$SampleType)  # unique items
        if(length(uniq) > 1 & ("B" %in% uniq)){
                stop("Your samples are a mix of children and adults. \n 
                       Split your samples in batches and rerun. \n
		       Look at Notes within this script for ideas on how to split")
        }
}

# ---- 2. Read data into RGChannelSets
message("\n")
message("Building RGChannelSet from array data...\n")

RGSet <- suppressWarnings(read.metharray.exp(targets=arrays, extended=FALSE))

# message("print size of session: \n")
# print(mem_used())

# ---- 4. Save the RGSet to disk
message("Saving the RGSet to disk...\n")

saveRDS(RGSet, file=output.rgset)

message("Building RGSet Done!")

