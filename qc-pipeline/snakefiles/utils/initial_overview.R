
debug <- T

if (debug) {
  
  args <- c(
    ""
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

id_file <- "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz"
id_data  <- read.table(
  file = id_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

fam_list <- lapply(args, function(file) {
  read.table(file, header = FALSE, stringsAsFactors = FALSE, col.names = c("fid", "iid", "pat", "mat", "sex", "phe"))
})

merged_fam <- do.call(rbind, fam_list)

merged_with_ids <- merged_fam %>% left_join(
    id_data %>% 
      select(
        iid = sentrix_id, ind_id = id, role
      ),
    by = "iid"
  )
