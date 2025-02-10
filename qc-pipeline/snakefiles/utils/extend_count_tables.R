debug <- F

if (debug) {
  
  args <- c(
    "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz",
    "",
    "/mnt/work/oystein/tmp/sample_counts",
    "mod1",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp001/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp002/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp003/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp007/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp008/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp009/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp010/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp011/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp012/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp014/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp015a/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp015b/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp016a/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp016b/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017a/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017b/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017c/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017d/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017e/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp017f/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp018a/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp018b/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp018c/initial_merge.fam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.07.01/mod1-data-preparation/snp018de/initial_merge.fam"
  )
  
} else {
  args <- commandArgs(TRUE)
}

library(janitor)
library(glue)
library(tidyr)
library(dplyr)
library(stringr)


id_file <- args[1]
existing_table <- args[2]
new_table <- args[3]
module <- args[4]
plink_files <- args[-c(1,2,3,4)]

get_batch_name <- function(file_path) {
  pattern <- "(snp\\d+[a-zA-Z]*)"
  batch <- str_extract(string = file_path, pattern = pattern)
  return(batch)
}

count_lines <- function(plink_file) {
  line_count <- as.integer(system(paste("wc -l <", plink_file), intern = TRUE))
  return(line_count)
}

count_list <- lapply(plink_files, count_lines)
batch_list <- lapply(plink_files, get_batch_name)

count_df <- data.frame(t(unlist(count_list)))
colnames(count_df) <- unlist(batch_list)
rownames(count_df) <- module

if(module != "mod1"){
    existing_df = read.table(existing_table)
    output_df <- rbind(existing_df, count_df)
} else {
    output_df <- count_df
}

write.table(output_df, new_table)


