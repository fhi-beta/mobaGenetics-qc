
##
#
# This script does a simple fam file reconstruction
#
##

set.seed(11111)


# Command line arguments
debug <- T
if (debug) {
  
  args <- c(
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/mod8_pedigree_ibd_estimate.kin0", 
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/mod8_check_sex.sexcheck",
    "/mnt/archive/snpQc/phenotypes/expected_relationship_24.04.12.gz",
    "/mnt/archive/snpQc/phenotypes/birth_year_24.04.12.gz",
    "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod7-post-imputation/all_samples/mod7_rename_missing_ids.psam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/mod8_psam_reconstruction.psam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/exclusion",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/mismatch_information.gz",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.09.04/mod8-release_annotation/mismatch_relationship.gz",
    "/mnt/work/oystein/tmp/fam_reconstruction_debug.md",
    "debug"
  )
  
} else {
 
args <- commandArgs(TRUE)

if (length(args) != 12) {
  
  stop(paste0("Nine arguments expected: genome file, sex check file, expected relationships file, birth year file, id file, current psam file, destination file, exclusion file, mismatch_information, mismatch_relationship, md file, title. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}
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


# Libraries

library(igraph)
library(dplyr)
library(ggplot2)
library(grid)

theme_set(theme_bw(base_size = 24))


# Load data

genomic_relatedness_table <- read.table(
  file = genome_file,
  header = T,
  stringsAsFactors = F
)
sex_check_data <- read.table(
  file = sex_check_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
expected_relationships_data <- read.table(
  file = expected_relationships_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
birth_year_data  <- read.table(
  file = birth_year_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

id_data  <- read.table(
  file = id_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

psam_data  <- read.table(
  file = psam_file,
  header = T,
  sep = "\t",
  col.names = c("IID", "SEX"),
  stringsAsFactors = F
)

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
  x = paste0("| Unknown | ", sum(mother_sex$SNPSEX == 0), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Male | ", sum(mother_sex$SNPSEX == 1), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Female | ", sum(mother_sex$SNPSEX == 2), " |\n"),
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
  x = paste0("| Unknown | ", sum(father_sex$SNPSEX == 0), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Male | ", sum(father_sex$SNPSEX == 1), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| Female | ", sum(father_sex$SNPSEX == 2), " |\n"),
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
    birth_year1 = ifelse(is.na(birth_year1), 1900, birth_year1),
    birth_year2 = ifelse(is.na(birth_year2), 1900, birth_year2)
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
  )

filtered_table <- subset(related_table, InfType =="PO")

parent_offspring_table <- data.frame(
    parent_sentrix_id = ifelse(filtered_table$birth_year1 < filtered_table$birth_year2, filtered_table$ID1,filtered_table$ID2),
    child_sentrix_id = ifelse(filtered_table$birth_year1 <filtered_table$birth_year2, filtered_table$ID2,filtered_table$ID1),
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
  psam_data,
  ids,
  id_type
) {

mother_offspring_detected <<- unique(parent_offspring_detected %>% filter(parent_sex == 2 & !is.na(parent_id) & !is.na(child_id) & age_difference> 12) %>% select(parent_id, child_id))
mother_offspring_detected_dupl <<- subset(mother_offspring_detected, duplicated(child_id))

father_offspring_detected <<- unique(parent_offspring_detected %>% filter(parent_sex == 1 & !is.na(parent_id) & !is.na(child_id) && age_difference > 12) %>% select(parent_id, child_id))
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

restored_psam_data <<- psam_data %>%
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


n_trios <<- nrow(restored_psam_data[!is.na(restored_psam_data$IID) & !is.na(restored_psam_data$MAT) & !is.na(restored_psam_data$PAT) & !duplicated(restored_psam_data),])
n_children <<- nrow(restored_psam_data[!is.na(restored_psam_data$IID) & (!is.na(restored_psam_data$MAT) | !is.na(restored_psam_data$PAT)) & !duplicated(restored_psam_data),])
n_fathers <<- length(unique(na.omit(restored_psam_data$PAT)))
n_mothers <<- length(unique(na.omit(restored_psam_data$MAT)))
n_mc <<- nrow(restored_psam_data[!is.na(restored_psam_data$IID) & !is.na(restored_psam_data$MAT) & !duplicated(restored_psam_data),])
n_fc <<- nrow(restored_psam_data[!is.na(restored_psam_data$IID) & !is.na(restored_psam_data$PAT) & !duplicated(restored_psam_data),])



write(
  x = paste("### ", tools::toTitleCase(id_type)),
  file = md_file,
  append = T
)
write(
  x = paste(n_distinct(ids), id_type,  "in total. Breakdown excluding multiple parents:\n - ", n_children, "children\n - ", n_mothers, "mothers\n - ", n_fathers, "fathers\n - ",  n_mc, "mother-child pairs\n - ", n_fc, "father-child pairs\n - ", n_trios, "trios\n"),
  file = md_file,
  append = T
)

write(
  x = paste("Multiple parents (at the", substring(id_type, 1, char(id_type)-1), "level):\n - ", nrow(mother_offspring_detected_dupl), "children with multiple mothers detected\n - ", nrow(father_offspring_detected_dupl), "children with multiple fathers detected\n - ", nrow(mother_offspring_expected_dupl), "children with multiple mothers in registry\n - ", nrow(father_offspring_expected_dupl), "children with multiple fathers in registry\n"),
  file = md_file,
  append = T
)


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
  x = paste(nrow(subset(psam_data, !(IID %in% id_data$sentrix_id))), "sentrix IDs missing from ID file"),
  file = md_file,
  append = T
)

check_expected_relationships(expected_relationships_data_ind, parent_offspring_table_ind, psam_data_ind, ind_ids, "individuals")
check_expected_relationships(expected_relationships_data_samples, parent_offspring_table_samples, psam_data, sample_ids, "samples")


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

 

