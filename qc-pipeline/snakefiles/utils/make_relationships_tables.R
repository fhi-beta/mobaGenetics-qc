library(janitor)
library(glue)
library(tidyr)
library(dplyr)

debug <- F

if (debug){
    args <- c("/mnt/archive2/moba_genotypes_resources/phenotypes/confirmed_relationships_25.01.30",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.01.30/mod8-release_annotation/mod8_pedigree_ibd_estimate.kin0",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_batch_table_batch",
    "/mnt/work/oystein/github/mobaGenetics-qc/qc-pipeline/snakefiles/parameters/batch_chip",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/ids_24.08.07.gz",
    "/mnt/archive3/phasing_test/expected_all_relations",
    "/mnt/archive3/phasing_test/expected_shapeit_fam",
    "/mnt/archive3/phasing_test/imputation_batches",
    "/mnt/archive3/phasing_test/expected_trios")
} else {
    args <- commandArgs(TRUE)
}

rel_file <- args[1]
gen_file <- args[2]
batch_file <- args[3]
chip_file <- args[4]
id_file <- args[5]
relations_file <- args[6]
shapeit_fam <- args[7]
imputation_batches_trunk <- args[8]
trios_file <- args[9]


ids <- read.table(id_file, header = T) %>% select(iid = sentrix_id, reg_id = id)
exp_rel <- read.table(rel_file, header = T, col.names=c("iid", "mat", "pat"))[,c(1,3,2)]
exp_rel <- subset(exp_rel, !is.na(iid) & (!is.na(pat) | !is.na(mat)))

gen <- read.table(gen_file, header = T)
po <- subset(gen, InfType == "PO")



# rel <- batches %>% select(iid) %>% left_join(exp_rel, by = "iid")

# rel <- subset(rel, !is.na(iid))
# rel <- subset(rel, !is.na(pat) | !is.na(mat))
batches <- read.table(batch_file, header = T)
batch_chip <- read.table(chip_file, header = T)

rel <- batches %>% select(iid) %>% left_join(exp_rel, by = "iid")


rel <- rel %>% left_join(po %>% select(iid = ID2, mat = ID1, mat1 = InfType), by = c("iid", "mat")) %>% left_join(po %>% select(iid = ID1, mat = ID2, mat2 = InfType), by = c("iid", "mat"))
rel$mat_confirmed <- ifelse(!is.na(rel$mat1) | !is.na(rel$mat2), 1, 0)

rel <- rel %>% left_join(po %>% select(iid = ID2, pat = ID1, pat1 = InfType), by = c("iid", "pat")) %>% left_join(po %>% select(iid = ID1, pat = ID2, pat2 = InfType), by = c("iid", "pat"))
rel$pat_confirmed <- ifelse(!is.na(rel$pat1) | !is.na(rel$pat2), 1, 0)

rel$pat <- ifelse(rel$pat_confirmed == 1, rel$pat, NA)
rel$mat <- ifelse(rel$mat_confirmed == 1, rel$mat, NA)

rel <- rel %>% select(iid, pat, mat)

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
rel <- rel %>% mutate(pat = ifelse(is.na(pat_batch), NA, pat)) %>% mutate(mat = ifelse(is.na(mat_batch), NA, mat))

rel$parents <- 0
rel$parents <- ifelse(!is.na(rel$pat), rel$parents +1, rel$parents)
rel$parents<- ifelse(!is.na(rel$mat), rel$parents +1, rel$parents)




rel_filtered <- rel %>%
  group_by(iid) %>%
  filter(shared_chips == max(shared_chips)) %>%
  ungroup() %>%
  as.data.frame()

rel_filtered <- rel_filtered %>%
  group_by(iid) %>%
  filter(parents_in_batch == max(parents_in_batch)) %>%
  ungroup() %>%
  as.data.frame()

rel_filtered <- rel_filtered %>%
  group_by(iid) %>%
  filter(parents == max(parents)) %>%
  ungroup() %>%
  as.data.frame()

rel_filtered <- rel_filtered %>% distinct(iid, .keep_all=TRUE)


shapeit_fam_df <- subset(rel_filtered, !is.na(pat) | !is.na(mat)) %>% select(iid, pat, mat)

write.table(
  x = rel_filtered,
  file =relations_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "NA"
)
trios <- subset(rel_filtered, !is.na(pat) & !is.na(mat))

write.table(
  x = trios,
  file =trios_file,
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "NA"
)


write.table(x = shapeit_fam_df, file = shapeit_fam, col.names = F, row.names = F, quote = F, sep = "\t")

batches <- unique(rel_filtered$iid_batch)

for (batch in batches) {
  samples_file <- paste0(imputation_batches_trunk, ".", batch)
  samples <- subset(rel_filtered, iid_batch == batch) %>% select(iid)
  samples <- unique(samples)
  write.table(samples, samples_file, row.names = F, col.names = F, quote=F)
}
