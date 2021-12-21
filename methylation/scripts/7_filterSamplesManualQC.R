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

message('Loading script dependencies...\n')


# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
})


message("printing arguments from snakemake...\n")
print(args)
message("\n")
# improve readability by unpacking args to variables:

input.methylSet <- args[1]
input.discarded_samples <- args[2]
input.filtered_samples = args[3]

output.methylSet = args[4]
output.removed_samples = args[5]

methylSet <-readRDS(input.methylSet)
discarded_samples <- readLines(input.discarded_samples)

# ---- 1. Filter methylset
message('Filtering out samples included in file:\n\t', input.discarded_samples, "\n")
all_samples = colnames(methylSet)
discarded_samples_logic = all_samples %in% discarded_samples

methylSetFiltered <- methylSet[,!discarded_samples_logic]

# ---- 2. Save methylset
message(paste0("Saving methylSet to ", output.methylSet, '...\n'))
saveRDS(methylSetFiltered, file=output.methylSet)

dropped_samples_df = read.csv(file = input.filtered_samples)
dropped_samples_df = rbind(dropped_samples_df, data.frame(sample=discarded_samples, reason=rep("manual_removed_bad_density", times = length(discarded_samples))))
write.csv(dropped_samples_df, file = output.removed_samples, row.names = F)

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("Manual filtering done!\n\n")

