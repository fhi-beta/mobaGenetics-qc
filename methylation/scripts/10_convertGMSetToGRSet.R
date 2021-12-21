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
    library(data.table, quietly=TRUE)
    library(pryr, quietly=TRUE)
})


message("Printing arguments from snakemake...\n")
print(args)
message('\n')

# Unpack args to improve readability:
input.methset = args[1]
output.ratioset = args[2]

# ---- 1. Read in data
message(paste0('\nReading in MethylSet from: ', input.methset, '...\n'))
methylSet <- readRDS(input.methset)


# ---- 2. Map MethylSet to GenomicMethylSet and save 
message('Mapping MethylSet to GenomicMethylSet...\n')

grSet <- mapToGenome(methylSet)
rm(methylSet)
invisible(gc())

# ---- 3. Convert GMSet to GRSet
message(paste0('\nConverting GenomicMethylSet to GenomicRatioSet...\n'))
grSet <- ratioConvert(grSet, what='both')

# ---- 4. Save to GenomicRatioSet
message(paste0('\nSaving GenomicRatioSet to: ', output.ratioset, '\n'))

saveRDS(grSet, output.ratioset)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("Converting to GRSet done!\n\n")
