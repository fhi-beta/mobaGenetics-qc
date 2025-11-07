library(janitor)
library(glue)
library(tidyr)
library(dplyr)

debug <- F

if (debug){
    args <- c("/mnt/archive2/moba_genotypes_resources/phenotypes/expected_relationship_24.04.12.gz",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_batch_table_batch",
    "/mnt/work/oystein/github/mobaGenetics-qc/qc-pipeline/snakefiles/parameters/batch_chip",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/ids_24.08.07.gz",
    "/mnt/archive3/phasing_test/expected_all_relations",
    "/mnt/archive3/phasing_test/expected_shapeit_fam")
} else {
    args <- commandArgs(TRUE)
}

rel_file <- args[1]
batch_file <- args[2]
chip_file <- args[3]
id_file <- args[4]
relations_file <- args[5]
shapeit_fam <- args[6]


ids <- read.table(id_file, header = T) %>% select(iid = sentrix_id, reg_id = id)
exp_rel <- read.table(rel_file, header = T, col.names=c("iid", "mat", "pat"))[,c(1,3,2)]
exp_rel <- subset(exp_rel, !is.na(iid) & (!is.na(pat) | !is.na(mat)))

write.table(x = exp_rel, file = shapeit_fam, col.names = F, row.names = F, quote = F, sep = "\t")

rel <- ids %>% select(iid) %>% left_join(exp_rel, by = "iid")

# rel <- subset(rel, !is.na(iid))
# rel <- subset(rel, !is.na(pat) | !is.na(mat))
batches <- read.table(batch_file, header = T)
batch_chip <- read.table(chip_file, header = T)

rel <- rel %>% left_join(batches %>% select(iid, iid_batch = batch), by = "iid")
rel <- rel %>% left_join(batches %>% select(pat = "iid", pat_batch = batch), by = "pat")
rel <- rel %>% left_join(batches %>% select(mat = "iid", mat_batch = batch), by = "mat")
rel <- rel %>% left_join(batch_chip %>% select(iid_batch = batch, iid_chip = chip), by ="iid_batch")
rel <- rel %>% left_join(batch_chip %>% select(pat_batch = batch, pat_chip = chip), by ="pat_batch")
rel <- rel %>% left_join(batch_chip %>% select(mat_batch = batch, mat_chip = chip), by ="mat_batch")
rel <- rel %>% left_join(ids %>% select(iid, iid_reg = reg_id), by = "iid")
rel <- rel %>% left_join(ids %>% select(pat = "iid", pat_reg = reg_id), by = "pat")
rel <- rel %>% left_join(ids %>% select(mat = "iid", mat_reg = reg_id), by = "mat")

rel$parents_in_batch <- 0
rel$parents_in_batch <- ifelse(!is.na(rel$pat_batch) & rel$iid_batch == rel$pat_batch, rel$parents_in_batch +1, rel$parents_in_batch)
rel$parents_in_batch <- ifelse(!is.na(rel$mat_batch) & rel$iid_batch == rel$mat_batch, rel$parents_in_batch +1, rel$parents_in_batch)
rel$shared_chips <- 0
rel$shared_chips <- ifelse(!is.na(rel$pat_chip) & rel$iid_chip == rel$pat_chip, rel$shared_chips +1, rel$shared_chips)
rel$shared_chips <- ifelse(!is.na(rel$mat_chip) & rel$iid_chip == rel$mat_chip, rel$shared_chips +1, rel$shared_chips)

rel <- subset(rel, !is.na(iid_batch))
write.table(
  x = rel,
  file =relations_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "NA"
)
