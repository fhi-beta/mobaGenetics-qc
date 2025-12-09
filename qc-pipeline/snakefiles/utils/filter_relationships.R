library(janitor)
library(glue)
library(tidyr)
library(dplyr)

debug <-F
if (debug) {
  args <- c("/mnt/archive2/moba_genotypes_resources/phenotypes/all_samples_confirmed_relationships_2025.09.25",
    "/mnt/archive3/phasing_test/filtered_shapeit_fam",
    "/mnt/archive3/phasing_test/filtered_relations",
    "/mnt/archive3/phasing_test/filtered_trios")
} else {
  args <- commandArgs(TRUE)
}

all_samples_file <- args[1]
shapeit_fam_file <- args[2]
relations_file <- args[3]
trios_file <- args[4]

rel <- read.table(all_samples_file, header = T)
rel$avg_parent_miss <- ifelse(is.na(rel$avg_parent_miss), 1, rel$avg_parent_miss)

filtered_rel <- rel %>%
  group_by(IID) %>%
  filter(shared_chips == max(shared_chips)) %>%
  ungroup() %>%
  as.data.frame()

filtered_rel <- filtered_rel %>%
  group_by(IID) %>%
  filter(parents_in_batch == max(parents_in_batch)) %>%
  ungroup() %>%
  as.data.frame()

filtered_rel <- filtered_rel %>%
  group_by(IID) %>%
  filter(expected_parents == max(expected_parents)) %>%
  ungroup() %>%
  as.data.frame()

filtered_rel <- filtered_rel %>%
  group_by(IID) %>%
  filter(avg_parent_miss == min(avg_parent_miss)) %>%
  ungroup() %>%
  as.data.frame()

filtered_rel <- filtered_rel %>% distinct(IID, .keep_all=TRUE)

filtered_rel <- subset(filtered_rel, !is.na(MAT) | !is.na(PAT))

shapeit_fam <- filtered_rel %>% select(IID, PAT, MAT)

trios <- subset(filtered_rel, !is.na(MAT) & !is.na(PAT)) %>% select(PAT, MAT, IID)



write.table(
  x = filtered_rel,
  file =relations_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "NA"
)

write.table(
  x = trios,
  file =trios_file,
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "NA"
)


write.table(x = shapeit_fam, file = shapeit_fam_file, col.names = F, row.names = F, quote = F, sep = "\t")

