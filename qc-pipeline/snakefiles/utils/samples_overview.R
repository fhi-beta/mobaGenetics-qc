debug <- F

if (debug) {
  
  args <- c(
    "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz",
    "/mnt/work/oystein/tmp/module_overview.md",
    "Before QC:",
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
md_file <- args[2]
title <- args[3]

write(
  x = title,
  file = md_file,
  append = T
)

id_data  <- read.table(
  file = id_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

write_counts <- function(fam_table){
  n_samples <- nrow(fam_table)
  n_individuals <- length(unique(fam_table$ind_id))
  n_samples_without_ind_id <- nrow(subset(fam_table, is.na(ind_id)))
  write(
    x = paste0(" - ", n_samples, " samples\n - ", n_individuals, " individuals\n - ", n_samples_without_ind_id, " samples without registry match\n\n"),
    file = md_file,
    append = T
)
}

join_and_write_counts <- function(file) {
  df <-  read.table(file, header = FALSE, stringsAsFactors = FALSE, col.names = c("fid", "iid", "pat", "mat", "sex", "phe")) %>% left_join(
    id_data %>% 
      select(
        iid = sentrix_id, ind_id = id, role
      ),
    by = "iid"
  )
  pattern <- "(snp\\d+)"
  batch <- str_extract(string = file, pattern = pattern)
  write(
    x = paste0("### ", batch),
    file = md_file,
    append = T
  )
  write_counts(df)
  return(df)
}

fam_list <- lapply(args[-c(1,2,3)], join_and_write_counts)

merged_fam <- do.call(rbind, fam_list)

 write(
    x = paste0("### Total:"),
    file = md_file,
    append = T
  )
write_counts(merged_fam)


