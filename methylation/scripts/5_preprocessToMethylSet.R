# ---- 0. Load dependencies

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

# activate renv if renv = TRUE
if(as.logical(args[1])){
        source("renv/activate.R")
}

args = args[-1]

message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
    library(pryr, quietly=TRUE)
})


message("printing arguments from snakemake...\n")
print(args)
message("\n")

# improve readability by unpacking args to variables:

input.rgset = args[1]

params.preprocess_methods = args[2]

output.methylsets = args[3]

# ---- 1. Read in rgSet
message(paste0('Reading input file from ', args[1], '...\n'))

rgSet <- readRDS(input.rgset)

for (i in seq_along(output.methylsets)) {

    # ---- 2. Preprocess the rgSet to one or more MethylSets using the selected preprocessing methods
    method <- paste0('preprocess', params.preprocess_methods)
    message(paste0("Preprocessing with ", method, " function...\n"))

    .preprocessFunction <- get(method)  # Gets the object named the same as the string passed to it

    if (method == 'preprocessFunnorm') {
        methylSet <- .preprocessFunction(rgSet, ratioConvert=FALSE)
    } else if (method == 'preprocessNoob') {
	methylSet <- .preprocessFunction(rgSet, dyeMethod="single")
    } else {
        methylSet <- .preprocessFunction(rgSet)
    }


    # ---- 3. Write output files
    message(paste0("Saving methylSet to ", output.methylsets[i], '...\n'))

    saveRDS(methylSet, file=output.methylsets[i])
    
    message("\n\n")
}


message("Preprocessing done!\n\n")
