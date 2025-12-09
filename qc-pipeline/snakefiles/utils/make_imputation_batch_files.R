library(dplyr)
debug <- F

if (debug) {
  args <- c("/mnt/archive3/phasing_test/phase_confirmed_relationships_chr3/mod6_batch_table",
  "/mnt/archive3/phasing_test/imputation_batches")
} else {
  args <- commandArgs(TRUE)
}

batch_table <- read.table(args[1], header = T)
imputation_batches_trunk <- args[2]
batches <- unique(batch_table$batch)
for (b in batches) {
  samples_file <- paste0(imputation_batches_trunk, ".", b)
  samples <- subset(batch_table, batch == b) %>% select(iid)
  samples <- unique(samples)
  write.table(samples, samples_file, row.names = F, col.names = F, quote=F)
}