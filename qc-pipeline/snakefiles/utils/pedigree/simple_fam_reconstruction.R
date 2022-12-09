
##
#
# This script does a simple fam file reconstruction
#
##


# Command line arguments
args <- commandArgs(TRUE)

if (length(args) != 8) {
  
  stop(paste0("Eight arguments expected: genome file, sex check file, registry file, current fam file, destination file, exclusion file, md file, title. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

genome_file <- args[1]

if (!file.exists(genome_file)) {
  
  stop("Genome file not found")
  
}

sex_check_file <- args[2]

if (!file.exists(sex_check_file)) {
  
  stop("sex check file not found")
  
}

registry_file <- args[3]

if (!file.exists(registry_file)) {
  
  stop("Registry file not found")
  
}

fam_file <- args[4]

if (!file.exists(fam_file)) {
  
  stop("Fam file not found")
  
}

destination_file <- args[5]

exclusion_file <- args[6]

md_file <- args[7]

title <- args[8]


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
registry_data  <- read.table(
  file = registry_file,
  header = T,
  sep = ",",
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
  family = paste0("family_", nuclear_families$membership),
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


# Map parents

mother_ids <- registry_data$SENTRIX_ID[registry_data$ROLE == "Mother"]
child_mother_table_1 <- related_table[related_table$InfType == "PO" & related_table$ID2 %in% mother_ids, c("ID1", "ID2")]
child_mother_table_2 <- related_table[related_table$InfType == "PO" & related_table$ID1 %in% mother_ids, c("ID1", "ID2")]
child_mother_table <- data.frame(
  id = c(child_mother_table_1[, 1], child_mother_table_1[, 2]),
  mother_id = c(child_mother_table_1[, 2], child_mother_table_1[, 1])
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

father_ids <- registry_data$SENTRIX_ID[registry_data$ROLE == "Father"]
child_father_table_1 <- related_table[related_table$InfType == "PO" & related_table$ID2 %in% father_ids, c("ID1", "ID2")]
child_father_table_2 <- related_table[related_table$InfType == "PO" & related_table$ID1 %in% father_ids, c("ID1", "ID2")]
child_father_table <- data.frame(
  id = c(child_father_table_1[, 1], child_father_table_1[, 2]),
  father_id = c(child_father_table_1[, 2], child_father_table_1[, 1])
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
  append = F,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)


# Write a list of samples to be excluded due to mismatch between the genetic relationship and registry information
child_ids <- registry_data$SENTRIX_ID[registry_data$ROLE == "Child"]
conflicts <- related_table[related_table$InfType == "PO" & related_table$ID1 %in% child_ids & related_table$ID2 %in% child_ids, c("ID1", "ID2")]
to_remove_ids <- c(conflicts$ID1, conflicts$ID2)
to_remove_fam <- updated_fam_data[updated_fam_data$id %in% to_remove_ids, c("family", "id")]

print(paste0("To remove: ", length(to_remove_ids)))
print(paste0("To remove: ", nrow(to_remove_fam)))

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
  x = paste0("- Number of children with parent-offspring relationship: ", nrow(to_remove_fam)),
  file = md_file,
  append = T
)


