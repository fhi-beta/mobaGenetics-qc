message("Load .Rprofile for set up... \n")
source(".Rprofile")

message("Loading script dependencies...\n")

# Suppres package load messages
suppressMessages({
        library(minfi, quietly=TRUE)
})

# Parse arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

message("Printing arguments from snakemake...\n")
print(args)
message('\n')

# Unpack args to improve readability:
input.rgset = args[1]
params.cellType = args[2]
params.preprocess = args[3]
output.cellCounts = args[4]

# ---- 1.Read in RGSet
message("\nReading in RGSet...\n")

RGSet <- readRDS(input.rgset)

# ---- 2. estimate cell counts

cellCounts = suppressWarnings(estimateCellCounts(RGSet, compositeCellType = params.cellType, 
						 processMethod = paste0("preprocess", params.preprocess), dyeMethod="single"))

write.csv(cellCounts, file = output.cellCounts, row.names = TRUE)

message(paste0("Cell counts estimated and stored: ", output.cellCounts, " \n\n"))

message('Cell count estimation done!')
