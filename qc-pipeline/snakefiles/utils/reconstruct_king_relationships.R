# Script for generating genetically confirmed relationships from KING output

library(igraph)
library(dplyr)
library(ggplot2)
library(grid) 
 
 args <- c(
    "/mnt/archive3/snpQc/pipeOut_dev/2025.01.30/mod8-release_annotation/mod8_pedigree_ibd_estimate.kin0", 
    "/mnt/archive3/snpQc/pipeOut_dev/2025.01.30/mod8-release_annotation/mod8_check_sex_common_filter_X_ycounts.sexcheck",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/birth_year_24.04.12.gz",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_batch_table_batch",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/confirmed_relationships_25.01.30")


genome_file <- args[1]

if (!file.exists(genome_file)) {
  
  stop("Genome file not found")
  
}

sex_check_file <- args[2]

if (!file.exists(sex_check_file)) {
  
  stop("sex check file not found")
  
}

birth_year_file <- args[3]

if (!file.exists(birth_year_file)) {
  
  stop("Birth year file not found")
  
}
batch_file <- args[4]
if (!file.exists(batch_file)) {
  
  stop("Batch file not found")
  
}

confirmed_relations_file <- args[5]

genomic_relatedness_table <- read.table(
  file = genome_file,
  header = T,
  stringsAsFactors = F
)
sex <- read.table(
  file = sex_check_file,
  header = F,
  col.names = c("IID", "PEDSEX", "SNPSEX", "STATUS", "F", "YCOUNT"),
  stringsAsFactors = F
)

sex$SNPSEX[is.na(sex$SNPSEX)] <- 0
sex$PEDSEX[is.na(sex$PEDSEX)] <- 0
sex$sex <- ifelse(sex$SNPSEX == "0", sex$PEDSEX, sex$SNPSEX)
batches <- read.table(batch_file, header = T)


expected_relationships_data <- read.table(
  file = expected_relationships_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

print("Read birth year file")
birth_year_data  <- read.table(
  file = birth_year_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

print("Read id file")
id_data  <- read.table(
  file = id_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)


related_table <- genomic_relatedness_table %>% 
  select(
    ID1, ID2, InfType
  ) %>% 
  filter(
    InfType == "PO"
  )

related_table <- related_table %>% 
  left_join(
    birth_year_data %>% 
      select(
        ID1 = sentrix_id,
        birth_year1 = birth_year
      ),
    by = "ID1"
  ) %>% 
  left_join(
    birth_year_data %>% 
      select(
        ID2 = sentrix_id,
        birth_year2 = birth_year
      ),
    by = "ID2"
  ) %>% 
  mutate(
    birth_year1 = ifelse(is.na(birth_year1), 1900, birth_year1),
    birth_year2 = ifelse(is.na(birth_year2), 1900, birth_year2)
  ) %>% left_join(sex %>% select(ID1 = IID, sex1 = sex), by = "ID1") %>% 
  left_join(sex %>% select(ID2 = IID, sex2 = sex), by = "ID2")

parent_offspring_table <- data.frame(
    parent_sentrix_id = as.character(ifelse(related_table$birth_year1 < related_table$birth_year2, related_table$ID1,related_table$ID2)),
    child_sentrix_id = as.character(ifelse(related_table$birth_year1 <related_table$birth_year2, related_table$ID2,related_table$ID1)),
    parent_birth_year = ifelse(related_table$birth_year1 < related_table$birth_year2,related_table$birth_year1,related_table$birth_year2),
    child_birth_year = ifelse(related_table$birth_year1 < related_table$birth_year2, related_table$birth_year2,related_table$birth_year1),
    age_difference = abs(related_table$birth_year1 - related_table$birth_year2),
    parent_sex = ifelse(related_table$birth_year1 < related_table$birth_year2, related_table$sex1,related_table$sex2),
    child_sex = ifelse(related_table$birth_year1 < related_table$birth_year2, related_table$sex2,related_table$sex1)
)

parent_offspring_table <- subset(parent_offspring_table, age_difference > 12)

parent_offspring_table_samples <- parent_offspring_table %>% rename(parent_id = parent_sentrix_id, child_id = child_sentrix_id)

father_offspring <- subset(parent_offspring_table_samples, parent_sex == 1)
mother_offspring <- subset(parent_offspring_table_samples, parent_sex == 2)

new_rel <- batches %>% select(child_sentrix_id = iid) %>% left_join(mother_offspring %>% select(child_sentrix_id = child_id, mother_sentrix_id = parent_id), by = "child_sentrix_id") %>% left_join(father_offspring %>% select(child_sentrix_id = child_id, father_sentrix_id = parent_id), by = "child_sentrix_id")
new_rel <- subset(new_rel, !is.na(mother_sentrix_id) | !is.na(father_sentrix_id))
write.table(new_rel, file = confirmed_relations_file, quote = F, row.names = F, col.names= T, sep = "\t")