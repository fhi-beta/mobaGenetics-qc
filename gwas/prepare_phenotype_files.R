##
#
# This script prepares phenotype files for Regenie
#
##

# Seed for random choice of sample

set.seed(05072023)

# Libraries

library(conflicted)
library(janitor)
library(glue)
library(dplyr)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Command line input

args <- commandArgs(TRUE)

if (length(args) != 3) {
  
  stop(paste0("Three command line arguments expected. ", length(args), " found."))
  
}

raw_pheno_file <- args[1]

if (!file.exists(raw_pheno_file)) {
  
  stop(paste0("Phenotype file ", raw_pheno_file, " not found."))
  
}

id_mapping_file <- args[2]

if (!file.exists(id_mapping_file)) {

  stop(paste0("Ids mapping file ", id_mapping_file, " not found."))

}

gwas_pheno_folder <- args[3]

if (!dir.exists(gwas_pheno_folder)) {
  
  stop(paste0("GWAS phenotypes folder ", gwas_pheno_folder, " not found."))
  
}

psam_file <- args[4] # "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod7-post-imputation/all_samples/mod7_rename_missing_ids.psam"

if (!file.exists(psam_file)) {
  
  stop(paste0("Fam file ", fam_file, " not found."))
  
}

mother_father_child_ids_file <- args[5] # "/mnt/archive/snpQc/phenotypes/expected_relationship_24.04.12.gz"

if (!file.exists(psam_fmother_father_child_ids_fileile)) {
  
  stop(paste0("Fam file ", mother_father_child_ids_file, " not found."))
  
}


# Load data

print(paste0(Sys.time(), " - Loading data"))

raw_phenotypes <- read.table(
  file = raw_pheno_file,
  header = T,
  sep = ";"
) %>%
    clean_names()

ids_mapping <- read.table(
  file = id_mapping_file,
  header = T,
  sep = "\t"
) %>%
    clean_names()

psam <- read.table(
  file = psam_file,
  header = F,
  sep = "\t"
) %>%
    clean_names()

names(psam) <- c("sentrix_id", "sex")

mother_father_child_ids <- read.table(
  file = mother_father_child_ids_file,
  header = T,
  sep = "\t"
) %>%
  clean_names()


# merge sentrix ids

phenotypes <- raw_phenotypes %>% 
  mutate(
    child_id = paste(preg_id_hdgb, barn_nr, sep = "_")
  ) %>% 
  left_join(
    ids_mapping %>% 
      filter(
        sentrix_id %in% psam$sentrix_id & role == "child"
      ) %>% 
      select(
        child_id = id,
        child_sentrix_id = sentrix_id
      ),
    by = "child_id",
    multiple = "any"
  ) %>% 
  left_join(
    mother_father_child_ids %>% 
      filter(
        !is.na(child_sentrix_id)
      ),
    by = "child_sentrix_id",
    multiple = "any"
  ) %>% 
  mutate(
    mother_sentrix_id = ifelse(mother_sentrix_id %in% psam$sentrix_id, mother_sentrix_id, NA),
    father_sentrix_id = ifelse(father_sentrix_id %in% psam$sentrix_id, father_sentrix_id, NA)
  )


# Set up data frames for GWAS

print(paste0(Sys.time(), " - Setting up for gwas"))

pheno_table_gwas_child <- phenotypes %>% 
  filter(
    !is.na(child_sentrix_id)
  ) %>% 
  select(
    FID = child_sentrix_id,
    IID = child_sentrix_id,
    where(is.numeric)
  )

pheno_table_gwas_mother <- phenotypes %>% 
  filter(
    !is.na(mother_sentrix_id)
  ) %>% 
  select(
    FID = mother_sentrix_id,
    IID = mother_sentrix_id,
    where(is.numeric)
  )

pheno_table_gwas_father <- phenotypes %>% 
  filter(
    !is.na(father_sentrix_id)
  ) %>% 
  select(
    FID = father_sentrix_id,
    IID = father_sentrix_id,
    where(is.numeric)
  )


# Write tables

print(paste0(Sys.time(), " - Export of tables to ", gwas_pheno_folder))

write.table(
  x = pheno_table_gwas_child,
  file = file.path(gwas_pheno_folder, "pheno_child"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_mother,
  file = file.path(gwas_pheno_folder, "pheno_mother"),
  row.names = F,
  col.names = T,
  quote = F
)

write.table(
  x = pheno_table_gwas_father,
  file = file.path(gwas_pheno_folder, "pheno_father"),
  row.names = F,
  col.names = T,
  quote = F
)



