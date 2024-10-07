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

debug <- F

if (!debug) {

args <- commandArgs(TRUE)

} else {

   args <- c(
   "/mnt/archive/snpQc/phenotypes/kvalitetssikring_vekt_hoyde_ny_24.09.06.csv",
   "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz",
   "/mnt/work/qc_genotypes/2024.07.01/moba_qc_gwas",
   "/mnt/archive/moba_genotypes_releases/2024.07.01/moba_genotypes_2024.07.01.psam",
   "/mnt/archive/snpQc/phenotypes/expected_relationship_24.04.12.gz",
   "/mnt/archive/moba_genotypes_releases/2024.07.01/batch/moba_genotypes_2024.07.01_batches"
   )

}

if (length(args) != 6) {
  
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

psam_file <- args[4]

if (!file.exists(psam_file)) {
  
  stop(paste0("Fam file ", fam_file, " not found."))
  
}

mother_father_child_ids_file <- args[5]

if (!file.exists(mother_father_child_ids_file)) {
  
  stop(paste0("Fam file ", mother_father_child_ids_file, " not found."))
  
}

batch_file <- args[6]

if (!file.exists(batch_file)) {

  stop(paste0("Batch file ", batch_file, " not found."))

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

names(psam) <- c("family_id", "sentrix_id", "maternal_id", "paternal_id", "sex")

mother_father_child_ids <- read.table(
  file = mother_father_child_ids_file,
  header = T,
  sep = "\t"
) %>%
  clean_names()

batch_table <- read.table(
  file = batch_file,
  header = T,
  sep = "\t"
) %>%
  clean_names()


# Exclusion criteria 
# Note - in an actual GWAS pipeline we would have a standalone pipeline for phenotypes taking care of outliers and harmonization. Since this is just for QC we do handling of the phenotypes directly in here and keep things to a minimum.

print(paste0(Sys.time(), " - Exclusion criteria"))

raw_phenotypes <- raw_phenotypes %>% 
  mutate(
    birth_weight = ifelse(vekt < 1000 | vekt > 10000, NA, vekt), 
    mother_height = ifelse(mor_hoyde < 100 | mor_hoyde > 250, NA, mor_hoyde), 
    father_height = ifelse(far_hoyde < 100 | far_hoyde > 250, NA, far_hoyde), 
    father_height_q1 = ifelse(far_hoyde_skjemafar < 100 | far_hoyde_skjemafar > 250, NA, far_hoyde_skjemafar), 
    father_height_q2 = ifelse(far_hoyde_farskjema2 < 100 | far_hoyde_farskjema2 > 250, NA, far_hoyde_farskjema2),
    father_height = ifelse(is.na(father_height) & !is.na(father_height_q1), father_height_q1, father_height), 
    father_height = ifelse(is.na(father_height) & !is.na(father_height_q2), father_height_q1, father_height)
  ) %>% 
  select(
    preg_id_hdgb, barn_nr, birth_weight, mother_height, father_height
  )


# merge sentrix ids

print(paste0(Sys.time(), " - Setting up identifiers"))

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

pheno_table_gwas_child <- psam %>% 
      select(
        family_id,
        child_sentrix_id = sentrix_id,
        sex
      ) %>% 
  inner_join(
    phenotypes %>% 
      filter(
        !is.na(child_sentrix_id)
      ),
    by = "child_sentrix_id",
    multiple = "any"
  ) %>% 
  select(
    FID = family_id,
    IID = child_sentrix_id,
    where(is.numeric)
  )

pheno_table_gwas_mother <- psam %>% 
      select(
        family_id,
        mother_sentrix_id = sentrix_id,
        sex
      ) %>% 
  inner_join(
    phenotypes %>% 
      filter(
        !is.na(mother_sentrix_id)
      ),
    by = "mother_sentrix_id",
    multiple = "any"
  ) %>% 
  mutate(
    sex = 2
  ) %>% 
  select(
    FID = family_id,
    IID = mother_sentrix_id,
    where(is.numeric)
  )

pheno_table_gwas_father <- psam %>% 
  select(
    family_id,
    father_sentrix_id = sentrix_id,
    sex
  ) %>% 
  inner_join(
    phenotypes %>% 
      filter(
        !is.na(father_sentrix_id)
      ),
    by = "father_sentrix_id",
    multiple = "any"
  ) %>% 
  mutate(
    sex = 1
  ) %>% 
  select(
    FID = family_id,
    IID = father_sentrix_id,
    where(is.numeric)
  )


# Add batch information as 1-hot encoding

print(paste0(Sys.time(), " - Annotating batches"))

batches <-  c("snp001", "snp002", "snp003", "snp007", "snp008", "snp009", "snp010", "snp011", "snp012", "snp014", "snp015a", "snp015b", "snp016a", "snp016b", "snp017a", "snp017b", "snp017c", "snp017d", "snp017e", "snp017f", "snp018a", "snp018b", "snp018c", "snp018de")

add_batch_columns <- function(df){
  
  df <- df %>% 
    left_join(
      batch_table %>% 
        rename(
          IID = iid
        ),
      by = "IID"
    )
  
  for (batch_name in batches) {
    
    df <- df %>% 
      mutate(
        !!sym(batch_name) := ifelse(batch == batch_name, 1, 0)
      )
    
  }
  
  df <- df %>% select(-batch)
  
  return(df)
  
}

pheno_table_gwas_child <- add_batch_columns(pheno_table_gwas_child)
pheno_table_gwas_mother <- add_batch_columns(pheno_table_gwas_mother)
pheno_table_gwas_father <- add_batch_columns(pheno_table_gwas_father)


# Write tables

print(paste0(Sys.time(), " - Exporting tables to ", gwas_pheno_folder))

if (length(unique(pheno_table_gwas_child$IID)) != nrow(pheno_table_gwas_child)) {
  
  stop("Non-unique sentrix id for children.")
  
}

write.table(
  x = pheno_table_gwas_child,
  file = file.path(gwas_pheno_folder, "pheno_child"),
  row.names = F,
  col.names = T,
  quote = F
)

if (length(unique(pheno_table_gwas_mother$IID)) != nrow(pheno_table_gwas_mother)) {
  
  stop("Non-unique sentrix id for mothers")
  
}

write.table(
  x = pheno_table_gwas_mother,
  file = file.path(gwas_pheno_folder, "pheno_mother"),
  row.names = F,
  col.names = T,
  quote = F
)

if (length(unique(pheno_table_gwas_father$IID)) != nrow(pheno_table_gwas_father)) {
  
  stop("Non-unique sentrix id for fathers")
  
}

write.table(
  x = pheno_table_gwas_father,
  file = file.path(gwas_pheno_folder, "pheno_father"),
  row.names = F,
  col.names = T,
  quote = F
)

pheno_table_gwas_parents <- rbind(pheno_table_gwas_mother, pheno_table_gwas_father)

if (length(unique(pheno_table_gwas_parents$IID)) != nrow(pheno_table_gwas_parents)) {
  
  stop("Non-unique sentrix id for parents")
  
}

write.table(
  x = pheno_table_gwas_parents,
  file = file.path(gwas_pheno_folder, "pheno_parent"),
  row.names = F,
  col.names = T,
  quote = F
)


# Write ids
# Note: this can be handled once and for all by the genotyping pipeline

ids_folder <- file.path(gwas_pheno_folder, "id")

print(paste0(Sys.time(), " - Exporting ids to ", ids_folder))

id_table_child <- pheno_table_gwas_child %>% select(FID, IID)

write.table(
  x = id_table_child,
  file = file.path(ids_folder, "children_id_plink"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = " "
)

id_table_mother <- pheno_table_gwas_mother %>% select(FID, IID)

write.table(
  x = id_table_mother,
  file = file.path(ids_folder, "mothers_id_plink"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = " "
)

id_table_father <- pheno_table_gwas_father %>% select(FID, IID)

write.table(
  x = id_table_father,
  file = file.path(ids_folder, "fathers_id_plink"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = " "
)

id_table_parents <- pheno_table_gwas_parents %>% select(FID, IID)

write.table(
  x = id_table_parents,
  file = file.path(ids_folder, "parents_id_plink"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = " "
)





