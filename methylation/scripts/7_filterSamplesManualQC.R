# ---- 0. Load dependencies

message('Loading script dependencies...\n')


# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality

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

message("Manual filtering done!\n\n")

