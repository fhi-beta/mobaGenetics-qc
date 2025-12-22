
##
#
# This script does a simple fam or psam file reconstruction. For fam input/otput, set plink_version to 1, for psam input/output, set plink_version to 2
#
##

set.seed(11111)


# Command line arguments
debug <- F
debug_plink_version <- 1
debug_batch <- "snp002"
if (debug) {
  if(debug_plink_version == 1){
    args <- c(
    paste0("/mnt/archive3/qc_genotypes/pipeOut_dev/2025.12.21/mod2-genetic-relationship/",debug_batch,"/pedigree_ibd_estimate.kin0"), 
    paste0("/mnt/archive3/qc_genotypes/pipeOut_dev/2025.12.21/mod2-genetic-relationship/", debug_batch, "/check_sex.sexcheck"),
    "/mnt/archive2/moba_genotypes_resources/phenotypes/expected_relationship_24.04.12.gz",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/birth_year_24.04.12.gz",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/ids_24.08.07.gz",
    paste0("/mnt/archive3/qc_genotypes/pipeOut_dev/2025.12.21/mod2-genetic-relationship/", debug_batch, "/callrate_permanent_removal.fam"),
    paste0("/mnt/work/oystein/tmp/fam_reconstruction/plink1/", debug_batch, "/fam_reconstruction.fam"),
    paste0("/mnt/work/oystein/tmp/fam_reconstruction/plink1/", debug_batch, "/exclusion"),
    paste0("/mnt/work/oystein/tmp/fam_reconstruction/plink1/", debug_batch, "/mismatch_information.gz"),
    paste0("/mnt/work/oystein/tmp/fam_reconstruction/plink1/", debug_batch, "/mismatch_relationship.gz"),
    paste0("/mnt/work/oystein/tmp/fam_reconstruction/plink1/", debug_batch, "/fam_reconstruction_debug.md"),
    "debug",
    "1"
  )
  } else if(debug_plink_version == 2){
    args <- c(
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_pedigree_ibd_estimate.kin0", 
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_import_ycounts.sexcheck",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/expected_relationship_24.04.12.gz",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/birth_year_24.04.12.gz",
    "/mnt/archive2/moba_genotypes_resources/phenotypes/ids_24.08.07.gz",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.01.30/mod8-release_annotation/mod8_split_rare_common.psam",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/mod8_psam_reconstruction.psam",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/exclusion",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/mismatch_information.gz",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/mismatch_relationship.gz",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/fam_reconstruction_debug.md",
    "debug",
    "2",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.01.30/mod8-release_annotation/mod8_batch_table_batch",
    "/mnt/work/oystein/github/mobaGenetics-qc/qc-pipeline/snakefiles/parameters/batch_chip",
    "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_best_snps.smiss",
    "/mnt/work/oystein/tmp/fam_reconstruction/plink2/debug/all_samples_relationships"
  )
  }
  
  
} else {
 
args <- commandArgs(TRUE)

# if (length(args) != 13) {
  
#   stop(paste0("13 arguments expected: IBD estimate file, sex check file, expected relationships file, birth year file, id file, current psam/fam file, destination file, exclusion file, mismatch_information, mismatch_relationship, md file, title, plink version (1 or 2) ", length(args), " found: ", paste(args, collapse = ", ")))

# }
}
genome_file <- args[1]

if (!file.exists(genome_file)) {
  
  stop("Genome file not found")
  
}

sex_check_file <- args[2]

if (!file.exists(sex_check_file)) {
  
  stop("sex check file not found")
  
}

expected_relationships_file <- args[3]

if (!file.exists(expected_relationships_file)) {
  
  stop("Expected relationship file not found")
  
}

birth_year_file <- args[4]

if (!file.exists(birth_year_file)) {
  
  stop("Birth year file not found")
  
}

id_file <- args[5]

if (!file.exists(id_file)) {
  
  stop("ID file not found")
  
}

psam_file <- args[6]

if (!file.exists(psam_file)) {
  
  stop("Psam file not found")
  
}

destination_file <- args[7]

exclusion_file <- args[8]

mismatch_information_file <- args[9]

mismatch_relationship_file <- args[10]

md_file <- args[11]

title <- args[12]
plink_version <- args[13]
if (plink_version == 2) {
  batch_file <- args[14]
  chip_file <- args[15]
  smiss_file <- args[16]
  all_samples_file <- args[17]
}


# Libraries

library(igraph)
library(dplyr)
library(ggplot2)
library(grid)

theme_set(theme_bw(base_size = 24))


# Load data
print("Read genome file")

genomic_relatedness_table <- read.table(
  file = genome_file,
  header = T,
  stringsAsFactors = F
)

print("Read sex check file")
if (plink_version == "1"){
sex_check_data <- read.table(
  file = sex_check_file,
  header = T,
  stringsAsFactors = F
)
} else {
  sex_check_data <- read.table(
  file = sex_check_file,
  header = F,
  col.names = c("IID", "PEDSEX", "SNPSEX", "STATUS", "F", "YCOUNT"),
  stringsAsFactors = F
)
}


sex_check_data$SNPSEX[is.na(sex_check_data$SNPSEX)] <- 0
sex_check_data$PEDSEX[is.na(sex_check_data$PEDSEX)] <- 0

print("Read relationship file")
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

if(plink_version == "1"){
print("Read fam file")
 psam_data_raw  <- read.table(
  file = psam_file,
  header = F,
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHE"),
  stringsAsFactors = F
)
psam_data <- psam_data_raw[,c("IID", "SEX")]

} else if(plink_version == "2"){
  print("Read psam file")
  psam_data  <- read.table(
  file = psam_file,
  header = F,
  col.names = c("IID", "SEX"),
  sep = "\t",
  stringsAsFactors = F
)
batches <- read.table(
  file = batch_file,
  header = T,
  stringsAsFactors = F
)
batch_chip <- read.table(
  file = chip_file,
  header = T,
  stringsAsFactors = F
)
smiss <- read.table(
  file = smiss_file,
  header = F, col.names = c("FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS"),
  stringsAsFactors = F
)
}

print("Done reading")
sample_ids <- psam_data[, 1]

genomic_relatedness_table$relationship <- factor(genomic_relatedness_table$InfType, levels = c("Dup/MZ", "PO", "FS", "2nd", "3rd", "4th", "UN"))
levels(genomic_relatedness_table$relationship) <- c("Duplicates or monozygotic twins", "Parent-offspring", "Full siblings", "2nd degree", "3rd degree", "4th degree", "Unrelated")

sentrix_not_in_mbr <- unique(sample_ids[!sample_ids %in% birth_year_data$sentrix_id])

if (length(sentrix_not_in_mbr) > 0) {
  
  mismatches_table <- data.frame(
    sentrix_id = sentrix_not_in_mbr,
    missing_birth_year = 1
  )
  
} else {
  
  mismatches_table <- data.frame(
    sentrix_id = character(),
    missing_birth_year = c()
  )
  
}

# Write docs

docs_dir <- dirname(md_file)
file_name <- basename(md_file)
docs_dir_name <- substr(file_name, 1, nchar(file_name) - 3)
docs_dir <- file.path(docs_dir, docs_dir_name)

if (!dir.exists(docs_dir)) {
  
  dir.create(docs_dir)
  
}

write(
  x = paste0("# ", title),
  file = md_file,
  append = F
)

write(
  x = paste0("- Number of samples in the genotyping data: ", nrow(psam_data), "."),
  file = md_file,
  append = T
)


write(
  x = "## Samples not in Medical Birth Regsitry",
  file = md_file,
  append = T
)
write(
  x = paste0(length(sentrix_not_in_mbr), " samples with missing birth year, assumed to be parent in the following."),
  file = md_file,
  append = T
)

write(
  x = "## Relationship inference",
  file = md_file,
  append = T
)

write(
  x = "| Relationship |   |",
  file = md_file,
  append = T
)
write(
  x = "| ------------ | - |",
  file = md_file,
  append = T
)

for (relationship in levels(genomic_relatedness_table$relationship)) {
  
  n <- sum(genomic_relatedness_table$relationship == relationship)
  write(
    x = paste0("| ", relationship, "| ", n, " |"),
    file = md_file,
    append = T
  )
}
write(
  x = "",
  file = md_file,
  append = T
)

ibd_plot <- ggplot() +
  geom_point(
    data = genomic_relatedness_table,
    mapping = aes(
      x = IBD1Seg,
      y = IBD2Seg,
      col = relationship
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "IBD1 segment",
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    name = "IBD2 segment",
    limits = c(0, 1)
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  ) + 
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1)
    )
  )

plot_name <- "ibd_plot.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(ibd_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

# Sex inference parents

mother_sex <- sex_check_data %>% 
  filter(
    IID %in% expected_relationships_data$mother_sentrix_id
  ) %>% 
  mutate(
    inferred_sex = factor(SNPSEX, levels = 0:2)
  )

  levels(mother_sex$inferred_sex) <- c("Unknown", "Male", "Female")

mother_exclude_sexcheck <- mother_sex$IID[mother_sex$SNPSEX == 1]

if (length(mother_exclude_sexcheck) > 0) {
  
  mismathes_mothers <- data.frame(
    sentrix_id = mother_exclude_sexcheck,
    mother_male = 1
  )
  
} else {
  
  mismathes_mothers <- data.frame(
    sentrix_id = character(),
    mother_male = c()
  )
  
}

mismatches_table <- mismatches_table %>% 
  full_join(
    mismathes_mothers,
    by = "sentrix_id"
  )

write(
  x = "## Mother sex check",
  file = md_file,
  append = T
)

write(
  x = "| Inferred sex |   |",
  file = md_file,
  append = T
)
write(
  x = "| ------------ | - |",
  file = md_file,
  append = T
)
write(
  x = paste0("| Unknown | ", nrow(subset(mother_sex, SNPSEX == 0)), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Male | ", nrow(subset(mother_sex, SNPSEX == 1)), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Female | ", nrow(subset(mother_sex, SNPSEX == 2)), " |\n"),
  file = md_file,
  append = T
)

mother_sex_plot <- ggplot() +
  geom_point(
    data = mother_sex,
    mapping = aes(
      x = F,
      y = YCOUNT,
      col = inferred_sex
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "F"
  ) +
  scale_y_continuous(
    name = "YCOUNT"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

plot_name <- "mother_sex_plot.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(mother_sex_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

father_sex <- sex_check_data %>% 
  filter(
    IID %in% expected_relationships_data$father_sentrix_id
  ) %>% 
  mutate(
    inferred_sex = factor(SNPSEX, levels = 0:2)
  )

levels(father_sex$inferred_sex) <- c("Unknown", "Male", "Female")

father_exclude_sexcheck <- father_sex$IID[father_sex$SNPSEX == 2]

if (length(father_exclude_sexcheck)) {
  
  father_sexcheck <- data.frame(
    sentrix_id = father_exclude_sexcheck,
    father_female = 1
  )
  
} else {
  
  father_sexcheck <- data.frame(
    sentrix_id = character(),
    father_female = c()
  )
  
}

mismatches_table <- mismatches_table %>% 
  full_join(
    father_sexcheck,
    by = "sentrix_id"
  )

write(
  x = "## Father sex check",
  file = md_file,
  append = T
)

write(
  x = "| Inferred sex |   |",
  file = md_file,
  append = T
)
write(
  x = "| ------------ | - |",
  file = md_file,
  append = T
)
write(
  x = paste0("| Unknown | ", nrow(subset(father_sex, SNPSEX == 0)), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Male | ", nrow(subset(father_sex, SNPSEX == 1)), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Female | ", nrow(subset(father_sex, SNPSEX == 2)), " |\n"),
  file = md_file,
  append = T
)

father_sex_plot <- ggplot() +
  geom_point(
    data = father_sex,
    mapping = aes(
      x = F,
      y = YCOUNT,
      col = inferred_sex
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "F"
  ) +
  scale_y_continuous(
    name = "YCOUNT"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

plot_name <- "father_sex_plot.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(father_sex_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

# Sex inference children

children_sex <- sex_check_data %>% 
  filter(
    ! IID %in% expected_relationships_data$mother_sentrix_id & ! IID %in% expected_relationships_data$father_sentrix_id
  ) %>% 
  mutate(
    inferred_sex = factor(SNPSEX, levels = 0:2)
  )

levels(children_sex$inferred_sex) <- c("Unknown", "Male", "Female")

write(
  x = "## Children sex check",
  file = md_file,
  append = T
)

write(
  x = "| Inferred sex |   |",
  file = md_file,
  append = T
)
write(
  x = "| ------------ | - |",
  file = md_file,
  append = T
)
write(
  x = paste0("| Unknown | ", sum(children_sex$SNPSEX == 0), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Male | ", sum(children_sex$SNPSEX == 1), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Female | ", sum(children_sex$SNPSEX == 2), " |\n"),
  file = md_file,
  append = T
)

children_sex_plot <- ggplot() +
  geom_point(
    data = children_sex,
    mapping = aes(
      x = F,
      y = YCOUNT,
      col = inferred_sex
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "F"
  ) +
  scale_y_continuous(
    name = "YCOUNT"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

plot_name <- "children_sex_plot.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(children_sex_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

# Make families

related_table <- genomic_relatedness_table %>% 
  select(
    ID1, ID2, InfType
  ) %>% 
  filter(
    InfType %in% c("Dup/MZ", "PO", "FS")
  )

relatedness_graph <- graph_from_data_frame(
  d = related_table,
  directed = F
)

nuclear_families <- components(relatedness_graph)

id_to_family <- data.frame(
  id = V(relatedness_graph)$name,
  family = paste0("kinship_cluster_", nuclear_families$membership),
  stringsAsFactors = F
)

id_to_family_singletons <- data.frame(
  id = sample_ids[! sample_ids %in% id_to_family$id],
  stringsAsFactors = F
) %>% 
  mutate(
    family = paste0("singleton_", row_number())
  )

id_to_family <- rbind(id_to_family, id_to_family_singletons)

# Add sex

id_to_family_sex <- id_to_family %>% 
  left_join(
    sex_check_data %>% 
      select(
        id = IID,
        sex = SNPSEX
      ),
    by = "id"
  ) %>% 
  mutate(
    sex = ifelse(id %in% expected_relationships_data$father_sentrix_id & sex == 0, 1, sex),
    sex = ifelse(id %in% expected_relationships_data$mother_sentrix_id & sex == 0, 2, sex)
  )

if (sum(is.na(id_to_family_sex$sex)) > 0) {
  
  stop("Some individuals miss sex inference")
  
}

# Merge with birth year and sex

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
    birth_year1 = ifelse(is.na(birth_year1), 0, birth_year1),
    birth_year2 = ifelse(is.na(birth_year2), 0, birth_year2)
  ) %>% 
  left_join(
    id_to_family_sex %>% 
      select(
        ID1 = id,
        sex1 = sex
      ),
    by = "ID1"
  ) %>% 
  left_join(
    id_to_family_sex %>% 
      select(
        ID2 = id,
        sex2 = sex
      ),
    by = "ID2"
  )%>% 
  left_join(
    id_data %>% 
      select(
        ID1 = sentrix_id,
        role1 = role
      ),
    by = "ID1"
  ) %>%
   left_join(
    id_data %>% 
      select(
        ID2 = sentrix_id,
        role2 = role
      ),
    by = "ID2"
  )

related_table$birth_year1 <- ifelse(related_table$birth_year1 == 0 & !is.na(related_table$role2) & related_table$role2 == "child", 1900, related_table$birth_year1)
related_table$birth_year1 <- ifelse(related_table$birth_year1 == 0 & !is.na(related_table$role2) & (related_table$role2 == "father" | related_table$role2 == "mother"), 2100, related_table$birth_year1)

related_table$birth_year2 <- ifelse(related_table$birth_year2 == 0 & !is.na(related_table$role1) & related_table$role1 == "child", 1900, related_table$birth_year2)
related_table$birth_year2 <- ifelse(related_table$birth_year2 == 0 & !is.na(related_table$role1) & (related_table$role1 == "father" | related_table$role1 == "mother"), 2100, related_table$birth_year2)


filtered_table <- subset(related_table, InfType =="PO")

parent_offspring_table <- data.frame(
    parent_sentrix_id = as.character(ifelse(filtered_table$birth_year1 < filtered_table$birth_year2, filtered_table$ID1,filtered_table$ID2)),
    child_sentrix_id = as.character(ifelse(filtered_table$birth_year1 <filtered_table$birth_year2, filtered_table$ID2,filtered_table$ID1)),
    parent_birth_year = ifelse(filtered_table$birth_year1 < filtered_table$birth_year2,filtered_table$birth_year1,filtered_table$birth_year2),
    child_birth_year = ifelse(filtered_table$birth_year1 < filtered_table$birth_year2, filtered_table$birth_year2,filtered_table$birth_year1),
    age_difference = abs(filtered_table$birth_year1 - filtered_table$birth_year2),
    parent_sex = ifelse(filtered_table$birth_year1 < filtered_table$birth_year2, filtered_table$sex1,filtered_table$sex2),
    child_sex = ifelse(filtered_table$birth_year1 < filtered_table$birth_year2, filtered_table$sex2,filtered_table$sex1)
)





parent_offspring_table_ind <- parent_offspring_table %>%
  left_join(id_data, by = c("child_sentrix_id" = "sentrix_id")) %>%
  rename(child_id = id) %>%
  left_join(id_data, by = c("parent_sentrix_id" = "sentrix_id")) %>%
  rename(parent_id = id) %>%
  select(parent_id, child_id, parent_birth_year, child_birth_year, age_difference, parent_sex, child_sex)

  parent_offspring_table_samples <- parent_offspring_table %>% rename(parent_id = parent_sentrix_id, child_id = child_sentrix_id)



expected_relationships_data_ind <- expected_relationships_data %>%
  left_join(id_data, by = c("child_sentrix_id" = "sentrix_id")) %>%
  rename(child_id = id) %>%
  left_join(id_data, by = c("mother_sentrix_id" = "sentrix_id")) %>%
  rename(mother_id = id) %>%
  left_join(id_data, by = c("father_sentrix_id" = "sentrix_id")) %>%
  rename(father_id = id) %>%
  select(child_id, mother_id, father_id)

expected_relationships_data_samples <- expected_relationships_data %>% rename(child_id = child_sentrix_id, mother_id = mother_sentrix_id, father_id = father_sentrix_id)

ind_ids <- id_data %>%
   filter(sentrix_id %in% sample_ids) %>%
   pull(id)
ind_ids <- unique(ind_ids)

psam_data_ind <- data.frame(IID = ind_ids)


if ( (nrow(parent_offspring_table_ind) != nrow(parent_offspring_table_samples)) |  (nrow(parent_offspring_table_ind) != nrow(parent_offspring_table))) {
  
  stop("Error when generating parent_offspring tables ")
  
}

if ( (nrow(expected_relationships_data_ind) != nrow(expected_relationships_data_samples)) |  (nrow(expected_relationships_data_ind) != nrow(expected_relationships_data))) {
  
  stop("Error when generating expected_relationships tables ")
  
}

expected_relationships_data_samples <- expected_relationships_data_samples[!duplicated(expected_relationships_data_samples),]
expected_relationships_data_ind <- expected_relationships_data_ind[!duplicated(expected_relationships_data_ind),]
parent_offspring_table_samples <- parent_offspring_table_samples[!duplicated(parent_offspring_table_samples),]
parent_offspring_table_ind <- parent_offspring_table_ind[!duplicated(parent_offspring_table_ind),]

found_in <- function(
  main_data,
  reference,
  ids
) {

  result <- main_data %>% 
    filter(
      !is.na(child_id) & 
      !is.na(parent_id) & 
      child_id %in% ids & 
      parent_id %in% ids
    ) %>%
    select(
      parent_id,
      child_id
    ) %>%
    mutate(
      found = ifelse(
        paste(parent_id, child_id) %in% 
          paste(
            reference$parent_id,
            reference$child_id
          ),
        TRUE,
        FALSE
      )
    )
  
  return(result)
}

write_relationship_docs <- function(
  total,
  n_found,
  n_not_found,
  header,
  body
){
  
write(
  x = paste(total, header),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", n_found, " (", round(100 * n_found/total, digits = 2), "%) ", body),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", n_not_found, " (", round(100 * n_not_found/total, digits = 2), "%) not ", body),
  file = md_file,
  append = T
)
write(
  x = "\n",
  file = md_file,
  append = T
)
}







check_expected_relationships <- function(
  parent_offspring_expected,
  parent_offspring_detected,
  orig_psam_data,
  ids,
  id_type
) {



mother_offspring_detected <<- unique(parent_offspring_detected %>% filter(parent_sex == 2 & !is.na(parent_id) & !is.na(child_id) & age_difference> 12) %>% select(parent_id, child_id))
mother_offspring_detected_dupl <<- subset(mother_offspring_detected, duplicated(child_id))
father_offspring_detected <<- unique(parent_offspring_detected %>% filter(parent_sex == 1 & !is.na(parent_id) & !is.na(child_id) & age_difference > 12) %>% select(parent_id, child_id))
father_offspring_detected_dupl <<- subset(father_offspring_detected, duplicated(child_id))



mother_offspring_expected <<- unique(parent_offspring_expected %>%
  select(parent_id = mother_id, child_id) %>%
  filter(!is.na(parent_id) & !is.na(child_id)))

mother_offspring_expected_dupl <<- subset(mother_offspring_expected, duplicated(child_id))

father_offspring_expected <<- unique(parent_offspring_expected %>%
  select(parent_id = father_id, child_id) %>%
  filter(!is.na(parent_id) & !is.na(child_id)))

father_offspring_expected_dupl <<- subset(father_offspring_expected, duplicated(child_id))

mother_offspring_expected_found <<- found_in(mother_offspring_expected, mother_offspring_detected, ids)
father_offspring_expected_found <<- found_in(father_offspring_expected, father_offspring_detected, ids)

mother_offspring_detected_found <<- found_in(mother_offspring_detected, mother_offspring_expected, ids)
father_offspring_detected_found <<- found_in(father_offspring_detected, father_offspring_expected, ids)

simplified_psam <<- orig_psam_data %>%
  left_join(
    father_offspring_detected %>%
      select(
        PAT = parent_id,
        IID = child_id
      ),
      by = "IID",
      multiple = "first"
  ) %>%
  left_join(
    mother_offspring_detected %>%
      select(
        MAT = parent_id,
        IID = child_id
      ),
      by = "IID",
      multiple = "first"
  )


n_trios <<- nrow(simplified_psam[!is.na(simplified_psam$IID) & !is.na(simplified_psam$MAT) & !is.na(simplified_psam$PAT) & !duplicated(simplified_psam),])
children <<- simplified_psam[!is.na(simplified_psam$IID) & (!is.na(simplified_psam$MAT) | !is.na(simplified_psam$PAT)) & !duplicated(simplified_psam),]$IID 
fathers <<- unique(na.omit(simplified_psam$PAT))
mothers <<- unique(na.omit(simplified_psam$MAT))
unrelated <<- subset(simplified_psam, !(IID %in% children) & !(IID %in% fathers) & !(IID %in% mothers))$IID 
n_children <<- length(children)
n_fathers <<- length(fathers)
n_mothers <<- length(mothers)
n_unrelated <<- length(unrelated)
n_mc <<- nrow(simplified_psam[!is.na(simplified_psam$IID) & !is.na(simplified_psam$MAT) & !duplicated(simplified_psam),])
n_fc <<- nrow(simplified_psam[!is.na(simplified_psam$IID) & !is.na(simplified_psam$PAT) & !duplicated(simplified_psam),])



write(
  x = paste("### ", tools::toTitleCase(id_type)),
  file = md_file,
  append = T
)
write(
  x = paste(n_distinct(ids), id_type,  "in total. Breakdown excluding multiple same-sex parents:\n - ", n_children, "children\n - ", n_mothers, "mothers\n - ", n_fathers, "fathers\n - ",  n_mc, "mother-child pairs\n - ", n_fc, "father-child pairs\n - ", n_trios, "mother-father-child trios\n - ", n_unrelated, "unrelated\n"),
  file = md_file,
  append = T
)
if(plink_version == "2"){
write(
  x = paste("Multiple same-sex parents (at the", substring(id_type, 1, nchar(id_type)-1), "level):\n - ", nrow(mother_offspring_detected_dupl), "children with more than one mother detected\n - ", nrow(father_offspring_detected_dupl), "children with more than one father detected\n - ", nrow(mother_offspring_expected_dupl), "children with more than one mother in registry\n - ", nrow(father_offspring_expected_dupl), "children with more than one father in registry\n"),
  file = md_file,
  append = T
)
}


write_relationship_docs(nrow(mother_offspring_expected_found), sum(mother_offspring_expected_found$found), sum(!mother_offspring_expected_found$found), "mother-child relationships expected.", "recovered by genetic relationships.")
write_relationship_docs(nrow(father_offspring_expected_found), sum(father_offspring_expected_found$found), sum(!father_offspring_expected_found$found), "father-child relationships expected.", "recovered by genetic relationships.")

write_relationship_docs(nrow(mother_offspring_detected_found), sum(mother_offspring_detected_found$found), sum(!mother_offspring_detected_found$found), "mother-child relationships detected.", "matched to registry.")
write_relationship_docs(nrow(father_offspring_detected_found), sum(father_offspring_detected_found$found), sum(!father_offspring_detected_found$found), "father-child relationships detected.", "matched to registry.")

}

write(
  x = paste("## Parental relationships"),
  file = md_file,
  append = T
)

write(
  x = paste(nrow(subset(psam_data, !(IID %in% id_data$sentrix_id))), "sentrix IDs missing from ID file. These are not counted as individuals."),
  file = md_file,
  append = T
)

check_expected_relationships(expected_relationships_data_ind, parent_offspring_table_ind, psam_data_ind, ind_ids, "individuals")
check_expected_relationships(expected_relationships_data_samples, parent_offspring_table_samples, psam_data, sample_ids, "samples")

if (plink_version== 1){
  new_psam <- simplified_psam %>% select(IID, PAT, MAT) %>% left_join(id_to_family_sex %>%
           select(
          FID = family,
            IID = id,
            SEX = sex
           ), by = "IID") %>% select(FID, IID, PAT, MAT, SEX)
} else {


rel <- psam_data %>% select(IID) %>% left_join(id_to_family_sex %>%
           select(
          FID = family,
            IID = id
           ), by = "IID")
rel <- rel %>%
  left_join(
    father_offspring_detected %>%
      select(
        PAT = parent_id,
        IID = child_id
      ),
      by = "IID"
  ) %>%
  left_join(
    mother_offspring_detected %>%
      select(
        MAT = parent_id,
        IID = child_id
      ),
      by = "IID"
  ) %>% select(IID, PAT, MAT, FID)



rel <- rel %>% left_join(id_to_family_sex %>% select(IID = id, SEX = sex), by = "IID")

rel <- rel %>% left_join(father_offspring_detected_found %>% select(IID = child_id, PAT = parent_id, expected_father = found), by = c("IID","PAT"))
rel <- rel %>% left_join(mother_offspring_detected_found %>% select(IID = child_id, MAT = parent_id, expected_mother = found), by = c("IID","MAT"))

rel <- rel %>% left_join(batches %>% select(IID = iid, iid_batch = batch), by = "IID")
rel <- rel %>% left_join(batches %>% select(PAT = "iid", pat_batch = batch), by = "PAT")
rel <- rel %>% left_join(batches %>% select(MAT = "iid", mat_batch = batch), by = "MAT")
rel <- rel %>% left_join(batch_chip %>% select(iid_batch = batch, iid_chip = chip), by ="iid_batch")
rel <- rel %>% left_join(batch_chip %>% select(pat_batch = batch, pat_chip = chip), by ="pat_batch")
rel <- rel %>% left_join(batch_chip %>% select(mat_batch = batch, mat_chip = chip), by ="mat_batch")
rel <- rel %>% left_join(id_data %>% select(IID = sentrix_id, iid_reg = id), by = "IID")
rel <- rel %>% left_join(id_data %>% select(PAT = sentrix_id, pat_reg = id), by = "PAT")
rel <- rel %>% left_join(id_data %>% select(MAT = sentrix_id, mat_reg = id), by = "MAT")

rel <- rel %>% left_join(smiss %>% select(IID, iid_miss = F_MISS), by = "IID")
rel <- rel %>% left_join(smiss %>% select(PAT = IID, pat_miss = F_MISS), by = "PAT")
rel <- rel %>% left_join(smiss %>% select(MAT = IID, mat_miss = F_MISS), by = "MAT")
rel$avg_parent_miss <- ifelse(is.na(rel$mat_miss) & is.na(rel$pat_miss), 
                              1,
                              ifelse(is.na(rel$mat_miss), 
                                     rel$pat_miss,
                                     ifelse(is.na(rel$pat_miss), 
                                            rel$mat_miss,
                                            (rel$pat_miss + rel$mat_miss) / 2)))


rel$parents_in_batch <- 0
rel$parents_in_batch <- ifelse(!is.na(rel$pat_batch) & rel$iid_batch == rel$pat_batch, rel$parents_in_batch +1, rel$parents_in_batch)
rel$parents_in_batch <- ifelse(!is.na(rel$mat_batch) & rel$iid_batch == rel$mat_batch, rel$parents_in_batch +1, rel$parents_in_batch)
rel$shared_chips <- 0
rel$shared_chips <- ifelse(!is.na(rel$pat_chip) & rel$iid_chip == rel$pat_chip, rel$shared_chips +1, rel$shared_chips)
rel$shared_chips <- ifelse(!is.na(rel$mat_chip) & rel$iid_chip == rel$mat_chip, rel$shared_chips +1, rel$shared_chips)

# rel <- subset(rel, !is.na(iid_batch))
rel <- rel %>% mutate(PAT = ifelse(is.na(pat_batch), NA, PAT)) %>% mutate(MAT = ifelse(is.na(mat_batch), NA, MAT))

rel$parents <- 0
rel$parents <- ifelse(!is.na(rel$PAT), rel$parents +1, rel$parents)
rel$parents <- ifelse(!is.na(rel$MAT), rel$parents +1, rel$parents)

rel$expected_parents <- 0
rel$expected_parents <- ifelse(!is.na(rel$expected_father) & rel$expected_father, rel$expected_parents +1, rel$expected_parents)
rel$expected_parents <- ifelse(!is.na(rel$expected_mother) & rel$expected_mother, rel$expected_parents +1, rel$expected_parents)


new_psam <- rel %>%
  group_by(IID) %>%
  filter(shared_chips == max(shared_chips)) %>%
  ungroup() %>%
  as.data.frame()

new_psam <- new_psam %>%
  group_by(IID) %>%
  filter(parents_in_batch == max(parents_in_batch)) %>%
  ungroup() %>%
  as.data.frame()

new_psam <- new_psam %>%
  group_by(IID) %>%
  filter(expected_parents == max(expected_parents)) %>%
  ungroup() %>%
  as.data.frame()

new_psam <- new_psam %>%
  group_by(IID) %>%
  filter(avg_parent_miss == min(avg_parent_miss)) %>%
  ungroup() %>%
  as.data.frame()

new_psam <- new_psam %>% distinct(IID, .keep_all=TRUE)
rel$avg_parent_miss <- ifelse(rel$avg_parent_miss == 1, NA, rel$avg_parent_miss)

new_psam <- new_psam %>% select(FID, IID, PAT, MAT, SEX)

}



# Conflicting relationships

conflicting_relationship_table <- parent_offspring_table %>% 
  filter(
age_difference <= 12
  ) %>% 
  select(
    child_sentrix_id, parent_sentrix_id, age_difference
  )
parent_offspring_detected_found <- rbind(mother_offspring_detected_found, father_offspring_detected_found)

conflicting_relationship_table <- conflicting_relationship_table %>% 
  full_join(
    mother_offspring_expected_found %>% 
      filter(
        !found
      ) %>% 
      select(
        child_sentrix_id = child_id,
        parent_sentrix_id = parent_id
      ) %>% 
      mutate(
        child_sentrix_id = as.character(child_sentrix_id),
        missing_mother_child_genetic_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  ) %>% 
  full_join(
    father_offspring_expected_found %>% 
      filter(
        !found
      ) %>% 
      select(
        child_sentrix_id = child_id,
        parent_sentrix_id = parent_id
      ) %>% 
      mutate(
        child_sentrix_id = as.character(child_sentrix_id),
        missing_father_child_genetic_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  )%>% 
  full_join(
    parent_offspring_detected_found %>% 
      filter(
        !found
      ) %>% 
      select(
        child_sentrix_id = child_id,
        parent_sentrix_id = parent_id
      ) %>% 
      mutate(
        child_sentrix_id = as.character(child_sentrix_id),
        missing_registry_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  )

mismatches_table <- mismatches_table %>%
   full_join(
    conflicting_relationship_table %>%
     mutate(conflicting_relationship_child = 1) %>%
     select(sentrix_id = child_sentrix_id, conflicting_relationship_child),
     by = "sentrix_id",
     multiple = "first"
    ) %>%
    full_join(
      conflicting_relationship_table %>%
     mutate(conflicting_relationship_parent = 1) %>%
     select(sentrix_id = parent_sentrix_id, conflicting_relationship_parent),
     by = "sentrix_id",
     multiple = "first"
    )


mismatches_table <- unique(mismatches_table)
 # Complete new psam
 
# restored_psam_data$SEX <- NULL
#  restored_psam_data <- restored_psam_data %>%
#       left_join(
#         id_to_family_sex %>%
#            select(
#             FID = family,
#             IID = id,
#             SEX = sex
#            ),
#         by = "IID"
#       )

# restored_psam_data <- restored_psam_data[, c("FID", setdiff(names(restored_psam_data), "FID"))]
if(plink_version == 2){
colnames(new_psam)[colnames(new_psam) == "FID"] <- "#FID"
}


if (!identical(psam_data$IID, new_psam$IID)){
  stop("Restored psam IIDs do not match original")
}

if(plink_version == "1"){
  new_psam$PHE <- -9
  write.table(
  x = new_psam,
  file = destination_file,
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "0"
)
} else if(plink_version == "2"){
  write.table(
  x = new_psam,
  file = destination_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F,
  na = "0"
)
write.table(
  x = rel,
  file = all_samples_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

}




write.table(
  x = mismatches_table,
  file = gzfile(mismatch_information_file),
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

write.table(
  x = conflicting_relationship_table,
  file = gzfile(mismatch_relationship_file),
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

# Write a list of samples to be excluded due to mismatch between the genetic relationship and registry information
child_ids <- expected_relationships_data$child_sentrix_id

mismatches_to_exclude <- mismatches_table$sentrix_id[!is.na(mismatches_table$mother_male)]
relationships_to_exclude1 <- conflicting_relationship_table %>% 
  filter(
    !is.na(age_difference) | !is.na(missing_mother_child_genetic_relationship)
  )
relationships_to_exclude2 <- conflicting_relationship_table %>% 
  filter(
    !is.na(missing_father_child_genetic_relationship) & child_sentrix_id %in% relationships_to_exclude1$child_sentrix_id
  )

to_remove_ids <- c(
  mismatches_to_exclude,
  relationships_to_exclude1$child_sentrix_id, relationships_to_exclude1$parent_sentrix_id,
  relationships_to_exclude2$child_sentrix_id, relationships_to_exclude2$parent_sentrix_id
  )

to_remove_ids <- unique(to_remove_ids)

to_remove_psam <- new_psam[new_psam$IID %in% to_remove_ids, 1:2]

write.table(
  x = to_remove_psam,
  file = exclusion_file,
  append = F,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)

write(
  x = "## Exclusion",
  file = md_file,
  append = T
)

write(
  x = paste0("- Number of samples excluded: ", nrow(to_remove_psam)),
  file = md_file,
  append = T
)
