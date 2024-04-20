
##
#
# This script does a simple fam file reconstruction
#
##

set.seed(11111)


# Command line arguments
debug <- F
if (debug) {
  
  args <- c(
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/pedigree_ibd_estimate.kin0", 
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/check_sex.sexcheck",
    "/mnt/archive/snpQc/phenotypes/expected_relationship_24.04.12.gz",
    "/mnt/archive/snpQc/phenotypes/birth_year_24.04.12.gz",
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/callrate_permanent_removal.fam",
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/fam_reconstruction.fam",
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/exclusion",
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/mismatch_information.gz",
    "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod2-genetic-relationship/snp012/mismatch_relationship.gz",
    "/mnt/work/marc/tmp/fam_reconstruction_debug.md",
    "debug"
  )
  
} else {
 
args <- commandArgs(TRUE)

if (length(args) != 11) {
  
  stop(paste0("Nine arguments expected: genome file, sex check file, expected relationships file, birth year file, current fam file, destination file, exclusion file, mismatch_information, mismatch_relationship, md file, title. ", length(args), " found: ", paste(args, collapse = ", ")))
  
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

fam_file <- args[5]

if (!file.exists(fam_file)) {
  
  stop("Fam file not found")
  
}

destination_file <- args[6]

exclusion_file <- args[7]

mismatch_information_file <- args[8]

mismatch_relationship_file <- args[9]

md_file <- args[10]

title <- args[11]


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
  stringsAsFactors = F
)
expected_relationships_data  <- read.table(
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
fam_data  <- read.table(
  file = fam_file,
  header = F,
  sep = " ",
  stringsAsFactors = F
)

sample_ids <- fam_data[, 2]

genomic_relatedness_table$relationship <- factor(genomic_relatedness_table$InfType, levels = c("Dup/MZ", "PO", "FS", "2nd", "3rd", "4th", "UN"))
levels(genomic_relatedness_table$relationship) <- c("Duplicates or monozygotic twins", "Parent-offspring", "Full siblings", "2nd degree", "3rd degree", "4th degree", "Unrelated")


# Exclude samples without birth year

sentrix_not_in_mbr <- unique(sample_ids[!sample_ids %in% birth_year_data$sentrix_id])

mismatches_table <- data.frame(
  sentrix_id = sentrix_not_in_mbr
)
mismatches_table$missing_birth_year <- 1


# Write docs

write(
  x = paste0("# ", title),
  file = md_file,
  append = F
)

docs_dir <- dirname(md_file)
file_name <- basename(md_file)
docs_dir_name <- substr(file_name, 1, nchar(file_name) - 3)
docs_dir <- file.path(docs_dir, docs_dir_name)

if (!dir.exists(docs_dir)) {
  
  dir.create(docs_dir)
  
}

write(
  x = "## Samples not in Medical Birth Regsitry",
  file = md_file,
  append = T
)
write(
  x = paste0(length(sentrix_not_in_mbr), " samples with missing birth year, will be assumed to be parent."),
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


# Sex assignment parents

mother_sex <- sex_check_data %>% 
  filter(
    IID %in% expected_relationships_data$mother_sentrix_id
  ) %>% 
  mutate(
    inferred_sex = factor(SNPSEX, levels = 0:2)
  )

levels(mother_sex$inferred_sex) <- c("Unknown", "Male", "Female")

mother_exclude_sexcheck <- mother_sex$IID[mother_sex$SNPSEX == 1]

mismathes_mothers <- data.frame(
  sentrix_id = mother_exclude_sexcheck
)
mismathes_mothers$mother_male <- 1

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

father_sexcheck <- data.frame(
  sentrix_id = father_exclude_sexcheck
)
father_sexcheck$father_female <- 1

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

# Conflicting relationships

conflicting_relationship_table <- related_table %>% 
  mutate(
    age_difference = abs(birth_year2 - birth_year1),
    child_sentrix_id = ifelse(birth_year1 < birth_year2, ID1, ID2),
    parent_sentrix_id = ifelse(birth_year1 < birth_year2, ID2, ID1)
  ) %>% 
  filter(
    InfType == "PO" & age_difference <= 12
  ) %>% 
  select(
    child_sentrix_id, parent_sentrix_id, age_difference
  )

child_mother_missing <- expected_relationships_data %>% 
  filter(
    !is.na(child_sentrix_id) & !is.na(mother_sentrix_id) & child_sentrix_id %in% sample_ids & mother_sentrix_id %in% sample_ids
  ) %>% 
  select(
    child_sentrix_id, mother_sentrix_id
  ) %>% 
  left_join(
    related_table %>% 
      filter(
        InfType == "PO"
      ) %>% 
      select(
        child_sentrix_id = ID1,
        mother_sentrix_id = ID2,
        genetic_relationship1 = InfType
      ),
    by = c("child_sentrix_id", "mother_sentrix_id")
  ) %>% 
  left_join(
    related_table %>% 
      filter(
        InfType == "PO"
      ) %>% 
      select(
        child_sentrix_id = ID2,
        mother_sentrix_id = ID1,
        genetic_relationship2 = InfType
      ),
    by = c("child_sentrix_id", "mother_sentrix_id")
  ) %>% 
  mutate(
    genetic_relationship = ifelse(!is.na(genetic_relationship1), "Found", "Missing"),
    genetic_relationship = ifelse(!is.na(genetic_relationship2), "Found", genetic_relationship)
  )

child_father_missing <- expected_relationships_data %>% 
  filter(
    !is.na(child_sentrix_id) & !is.na(father_sentrix_id) & child_sentrix_id %in% sample_ids & father_sentrix_id %in% sample_ids
  ) %>% 
  select(
    child_sentrix_id, father_sentrix_id
  ) %>% 
  left_join(
    related_table %>% 
      filter(
        InfType == "PO"
      ) %>% 
      select(
        child_sentrix_id = ID1,
        father_sentrix_id = ID2,
        genetic_relationship1 = InfType
      ),
    by = c("child_sentrix_id", "father_sentrix_id")
  ) %>% 
  left_join(
    related_table %>% 
      filter(
        InfType == "PO"
      ) %>% 
      select(
        child_sentrix_id = ID2,
        father_sentrix_id = ID1,
        genetic_relationship2 = InfType
      ),
    by = c("child_sentrix_id", "father_sentrix_id")
  ) %>% 
  mutate(
    genetic_relationship = ifelse(!is.na(genetic_relationship1), "Found", "Missing"),
    genetic_relationship = ifelse(!is.na(genetic_relationship2), "Found", genetic_relationship)
  )

write(
  x = "## Parental relationship",
  file = md_file,
  append = T
)
write(
  x = paste0(nrow(child_mother_missing), " mother-child relationships expected."),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(child_mother_missing$genetic_relationship == "Found"), " (", round(100 * sum(child_mother_missing$genetic_relationship == "Found")/nrow(child_mother_missing), digits = 2), "%) recovered by genetic relationships."),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(child_mother_missing$genetic_relationship == "Missing"), " (", round(100 * sum(child_mother_missing$genetic_relationship == "Missing")/nrow(child_mother_missing), digits = 2), "%) not recovered by genetic relationships."),
  file = md_file,
  append = T
)
write(
  x = paste0(nrow(child_father_missing), " father-child relationships expected."),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(child_father_missing$genetic_relationship == "Found"), " (", round(100 * sum(child_father_missing$genetic_relationship == "Found")/nrow(child_father_missing), digits = 2), "%) recovered by genetic relationships."),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(child_father_missing$genetic_relationship == "Missing"), " (", round(100 * sum(child_father_missing$genetic_relationship == "Missing")/nrow(child_father_missing), digits = 2), "%) not recovered by genetic relationships."),
  file = md_file,
  append = T
)

conflicting_relationship_table <- conflicting_relationship_table %>% 
  full_join(
    child_mother_missing %>% 
      filter(
        genetic_relationship == "Missing"
      ) %>% 
      select(
        child_sentrix_id,
        parent_sentrix_id = mother_sentrix_id
      ) %>% 
      mutate(
        missing_mother_child_genetic_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  ) %>% 
  full_join(
    child_father_missing %>% 
      filter(
        genetic_relationship == "Missing"
      ) %>% 
      select(
        child_sentrix_id,
        parent_sentrix_id = father_sentrix_id
      ) %>% 
      mutate(
        missing_father_child_genetic_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  )

genetic_relationship_match <- related_table %>% 
  filter(
    InfType == "PO"
  ) %>% 
  left_join(
    expected_relationships_data %>% 
      select(
        ID1 = child_sentrix_id,
        ID2 = mother_sentrix_id,
        father_sentrix_id
      ) %>% 
      filter(
        !is.na(ID1) & !is.na(ID2)
      ) %>% 
      mutate(
        genetic_relationship1 = "child-mother"
      ),
    by = c("ID1", "ID2"),
    multiple = "all"
  ) %>% 
  left_join(
    expected_relationships_data %>% 
      select(
        ID1 = mother_sentrix_id,
        ID2 = child_sentrix_id
      ) %>% 
      filter(
        !is.na(ID1) & !is.na(ID2)
      ) %>% 
      mutate(
        genetic_relationship2 = "mother-child"
      ),
    by = c("ID1", "ID2"),
    multiple = "all"
  ) %>% 
  left_join(
    expected_relationships_data %>% 
      select(
        ID1 = child_sentrix_id,
        ID2 = father_sentrix_id
      ) %>% 
      filter(
        !is.na(ID1) & !is.na(ID2)
      ) %>% 
      mutate(
        genetic_relationship3 = "child-father"
      ),
    by = c("ID1", "ID2"),
    multiple = "all"
  ) %>% 
  left_join(
    expected_relationships_data %>% 
      select(
        ID1 = father_sentrix_id,
        ID2 = child_sentrix_id
      ) %>% 
      filter(
        !is.na(ID1) & !is.na(ID2)
      ) %>% 
      mutate(
        genetic_relationship4 = "father-child"
      ),
    by = c("ID1", "ID2"),
    multiple = "all"
  ) %>% 
  mutate(
    mismatch = ifelse(is.na(genetic_relationship1) & is.na(genetic_relationship2) & is.na(genetic_relationship3) & is.na(genetic_relationship4), "Missing", "Found")
  )

write(
  x = paste0(nrow(genetic_relationship_match), " parent-offspring relationships detected"),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(genetic_relationship_match$mismatch == "Found"), " (", round(100 * sum(genetic_relationship_match$mismatch == "Found")/nrow(genetic_relationship_match), digits = 2), "%) match to registry."),
  file = md_file,
  append = T
)
write(
  x = paste0("- ", sum(genetic_relationship_match$mismatch == "Missing"), " (", round(100 * sum(genetic_relationship_match$mismatch == "Missing")/nrow(genetic_relationship_match), digits = 2), "%) do not match to registry."),
  file = md_file,
  append = T
)

missing_genetic_relationship <- genetic_relationship_match %>% 
  filter(
    mismatch == "Missing"
  ) %>% 
  mutate(
    child_sentrix_id = ifelse(birth_year1 < birth_year2, ID1, ID2),
    parent_sentrix_id = ifelse(birth_year1 < birth_year2, ID2, ID1)
  )

conflicting_relationship_table <- conflicting_relationship_table %>% 
  full_join(
    missing_genetic_relationship %>% 
      select(
        child_sentrix_id,
        parent_sentrix_id
      ) %>% 
      mutate(
        missing_registry_relationship = 1
      ),
    by = c("child_sentrix_id", "parent_sentrix_id")
  )

mismatches_table <- mismatches_table %>% 
  mutate(
    conflicting_relationship_child = ifelse(sentrix_id %in% conflicting_relationship_table$child_sentrix_id, 1, NA),
    conflicting_relationship_parent = ifelse(sentrix_id %in% conflicting_relationship_table$parent_sentrix_id, 1, NA)
  )


# Map parents

child_mother_table_1 <- related_table[related_table$InfType == "PO" & related_table$sex2 == 2 & related_table$birth_year2 > related_table$birth_year1 + 12, c("ID1", "ID2")]
child_mother_table_2 <- related_table[related_table$InfType == "PO" & related_table$sex1 == 2 & related_table$birth_year1 > related_table$birth_year2 + 12, c("ID1", "ID2")]
child_mother_table <- data.frame(
  id = c(child_mother_table_1$ID1, child_mother_table_2$ID2),
  mother_id = c(child_mother_table_1$ID2, child_mother_table_2$ID1)
) %>% 
  distinct() %>% 
  group_by(id) %>% 
  arrange(
    mother_id
  ) %>% 
  filter(
    row_number() == 1
  ) %>% 
  ungroup()

id_to_family_sex_mother <- id_to_family_sex %>% 
  left_join(
    child_mother_table,
    by = "id"
  ) %>% 
  mutate(
    mother_id = ifelse(is.na(mother_id), 0, mother_id)
  )

child_father_table_1 <- related_table[related_table$InfType == "PO" & related_table$sex2 == 1 & related_table$birth_year2 > related_table$birth_year1 + 12, c("ID1", "ID2")]
child_father_table_2 <- related_table[related_table$InfType == "PO" & related_table$sex1 == 1 & related_table$birth_year1 > related_table$birth_year2 + 12, c("ID1", "ID2")]
child_father_table <- data.frame(
  id = c(child_father_table_1$ID1, child_father_table_2$ID2),
  father_id = c(child_father_table_1$ID2, child_father_table_2$ID1)
) %>% 
  distinct() %>% 
  group_by(id) %>% 
  arrange(
    father_id
  ) %>% 
  filter(
    row_number() == 1
  ) %>% 
  ungroup()

id_to_family_sex_mother_father <- id_to_family_sex_mother %>% 
  left_join(
    child_father_table,
    by = "id"
  ) %>% 
  mutate(
    father_id = ifelse(is.na(father_id), 0, father_id)
  )


# build fam file

updated_fam_data <- id_to_family_sex_mother_father %>% 
  select(
    family,
    id,
    father_id,
    mother_id,
    sex
  ) %>% 
  mutate(
    pheno = -9
  )

if (length(unique(updated_fam_data$id)) != nrow(updated_fam_data)) {
  
  stop("Non-unique identifier introduced in fam file.")
  
}


# Save

write.table(
  x = updated_fam_data,
  file = destination_file,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)

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
to_remove_fam <- updated_fam_data[updated_fam_data$id %in% to_remove_ids, c("family", "id")]

write.table(
  x = to_remove_fam,
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
  x = paste0("- Number of samples excluded: ", nrow(to_remove_fam)),
  file = md_file,
  append = T
)


