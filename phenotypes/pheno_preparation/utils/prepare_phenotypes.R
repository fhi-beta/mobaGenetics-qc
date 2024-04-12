

# Libraries

library(conflicted)
library(foreign)
library(stringr)
library(glue)
library(tidyr)
library(dplyr)
library(janitor)

conflicts_prefer(dplyr::filter)


# Path to files

pheno_folder <- "/mnt/archive2/MomicsSource/snpArray/qcDevops78030/Nye_2024_04_05"
mfr_variables <- "2024_04_05_MFRVariablar.sav"
preg_id_child <- "2024_04_05_MobaGenetics_PREGID_Child.sav"
preg_id_father <- "2024_04_05_MobaGenetics_PREGID_Father.sav"
preg_id_mother <- "2024_04_05_MobaGenetics_PREGID_Mother.sav"
keys <- "2024_04_05_NOKLER_PDBHDGB.sav"
smoking <- "2024_04_05_SmokingStatus.sav"

birth_year_file <- "/mnt/archive/snpQc/phenotypes/birth_year_24.04.12.gz"
expected_relationship_file <- "/mnt/archive/snpQc/phenotypes/expected_relationship_24.04.12.gz"


# Load data

sav_mfr <- read.spss(
  file = file.path(pheno_folder, mfr_variables), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_mfr)) {
  
  if (is.character(sav_mfr[[colName]])) {
    
    sav_mfr[[colName]] <- str_trim(sav_mfr[[colName]])
    
  }
}

sav_preg_id_mother <- read.spss(
  file = file.path(pheno_folder, preg_id_mother), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_preg_id_mother)) {
  
  if (is.character(sav_preg_id_mother[[colName]])) {
    
    sav_preg_id_mother[[colName]] <- str_trim(sav_preg_id_mother[[colName]])
    
  }
}

sav_preg_id_father <- read.spss(
  file = file.path(pheno_folder, preg_id_father), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_preg_id_father)) {
  
  if (is.character(sav_preg_id_father[[colName]])) {
    
    sav_preg_id_father[[colName]] <- str_trim(sav_preg_id_father[[colName]])
    
  }
}

sav_preg_id_child <- read.spss(
  file = file.path(pheno_folder, preg_id_child), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_preg_id_child)) {
  
  if (is.character(sav_preg_id_child[[colName]])) {
    
    sav_preg_id_child[[colName]] <- str_trim(sav_preg_id_child[[colName]])
    
  }
}

sav_smoking <- read.spss(
  file = file.path(pheno_folder, smoking), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_smoking)) {
  
  if (is.character(sav_smoking[[colName]])) {
    
    sav_smoking[[colName]] <- str_trim(sav_smoking[[colName]])
    
  }
}

sav_keys <- read.spss(
  file = file.path(pheno_folder, keys), 
  use.value.labels = T, 
  to.data.frame = T, 
  stringsAsFactors = F
) %>% 
  clean_names()

for (colName in names(sav_keys)) {
  
  if (is.character(sav_keys[[colName]])) {
    
    sav_keys[[colName]] <- str_trim(sav_keys[[colName]])
    
  }
}


# Get a data frame of birth years

birth_year_child <- data.frame(
  id = sav_preg_id_child$sentrix_id,
  birth_year = sav_preg_id_child$faar
) %>% 
  filter(
    !is.na(birth_year)
  ) %>% 
  distinct()

birth_year_mother <- data.frame(
  id = sav_keys$m_id_hdgb,
  birth_year = sav_keys$mor_faar
) %>% 
  distinct() %>% 
  filter(
    !is.na(birth_year)
  ) %>% 
  left_join(
    sav_preg_id_mother %>% 
      select(
        id = m_id_hdgb,
        sentrix_id
      ),
    by = "id"
  )

birth_year_father <- data.frame(
  id = sav_keys$f_id_hdgb,
  birth_year = sav_keys$far_faar
) %>% 
  distinct() %>% 
  filter(
    !is.na(birth_year)
  ) %>% 
  left_join(
    sav_preg_id_father %>% 
      select(
        id = f_id_hdgb,
        sentrix_id
      ),
    by = "id"
  )

birth_year_table <- rbind(birth_year_child, birth_year_mother, birth_year_father) %>% 
  filter(
    !is.na(sentrix_id)
  ) %>% 
  select(
    sentrix_id, birth_year
  )

if (length(unique(birth_year_table$sentrix_id)) != nrow(birth_year_table)) {
  
  stop("duplicate ids in birth year table")
  
}

write.table(
  x = birth_year_table,
  file = gzfile(birth_year_file),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)


# Get expected relationships

expected_relationship <- sav_keys %>%
  left_join(
    sav_preg_id_mother %>% 
      select(
        m_id_hdgb,
        mother_sentrix_id = sentrix_id
      ),
    by = "m_id_hdgb",
    relationship = "many-to-many"
  ) %>% 
  left_join(
    sav_preg_id_father %>% 
      select(
        f_id_hdgb,
        father_sentrix_id = sentrix_id
      ),
    by = "f_id_hdgb",
    relationship = "many-to-many"
  ) %>% 
  left_join(
    sav_preg_id_child %>% 
      select(
        preg_id_hdgb,
        child_sentrix_id = sentrix_id
      ),
    by = "preg_id_hdgb",
    relationship = "many-to-many"
  ) %>% 
  filter(
    !is.na(mother_sentrix_id) | !is.na(father_sentrix_id) | !is.na(child_sentrix_id)
  ) %>% 
  select(
    child_sentrix_id, mother_sentrix_id, father_sentrix_id
  )

write.table(
  x = expected_relationship,
  file = gzfile(expected_relationship_file),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)



