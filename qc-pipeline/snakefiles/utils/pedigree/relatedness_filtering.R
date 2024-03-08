
##
#
# This script builds sets of unrelated participants.
# It is based on the `accumPIHAT` script written by Jonas on 2017.01.15.
#
##

set.seed(11111)


# Command line arguments
debug <- T

if (debug) {
  
  king_file <- "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod4-good_markers/snp012/founders/ibd_estimate.kin0"
  
  fam_file <- "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod4-good_markers/snp012/founders/ibd_estimate.fam"
  
  kinship_threshold <- 0.1
  
  accumulated_kinship_threshold <- 0.015
  
  exclusion_file <- "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod4-good_markers/snp012/founders/ibd_accum_pedigree_sample_exclusion"
  
  md_file <- "/mnt/work/marc/github/mobaGenetics-qc/qc-pipeline/docs/snp012/ibd_accum/founders/ibd_accum.md"
  
  title <- "Kinship filtering in snp012"
  
} else {

args <- commandArgs(TRUE)

if (length(args) != 7) {
  
  stop(paste0("Six arguments expected. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

king_file <- args[1]

if (!file.exists(king_file)) {
  
  stop("King file not found")
  
}

fam_file <- args[2]

if (!file.exists(fam_file)) {
  
  stop("Fam file not found")
  
}

kinship_threshold <- as.numeric(args[3])

if (is.na(kinship_threshold)) {
  
  stop(psate0("Input for `kinship_threshold` (", kinship_threshold, ") could not be parsed as a number."))
  
}

accumulated_kinship_threshold <- as.numeric(args[4])

if (is.na(accumulated_kinship_threshold)) {
  
  stop(psate0("Input for `accumulated_kinship_threshold` (", accumulated_kinship_threshold, ") could not be parsed as a number."))
  
}

exclusion_file <- args[5]

md_file <- args[6]

title <- args[7]

}


# Libraries

library(glue)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggside)
library(grid)

theme_set(theme_bw(base_size = 24))


# Load data

print(glue("{Sys.time()}    Relatedness analysis - Loading data"))

genomic_relatedness_table <- read.table(
  file = king_file,
  header = T,
  stringsAsFactors = F
)

fam_data  <- read.table(
  file = fam_file,
  header = F,
  sep = " ",
  stringsAsFactors = F
)

ids <- unique(fam_data[, 2])


# Write docs

print(glue("{Sys.time()}    Relatedness analysis - Writing documentation"))

write(
  x = paste0("# ", title),
  file = md_file,
  append = F
)

write(
  x = "Relatedness filtering, {length(ids)} individuals.",
  file = md_file,
  append = T
)

docs_dir <- dirname(md_file)
file_name <- basename(md_file)
docs_dir_name <- substr(file_name, 1, nchar(file_name) - 3)
docs_dir <- file.path(docs_dir, docs_dir_name)

if (!dir.exists(docs_dir)) {
  
  dir.create(docs_dir)
  
}

write(
  x = "## Relatedness",
  file = md_file,
  append = T
)

ibd_plot <- ggplot() +
  geom_point(
    data = genomic_relatedness_table,
    mapping = aes(
      x = IBS0,
      y = Kinship
    ),
    alpha = 0.2
  ) +
  geom_ysidedensity(
    data = genomic_relatedness_table,
    mapping = aes(
      x = after_stat(density),
      y = Kinship
    ),
    fill = "grey80"
  ) +
  geom_hline(
    yintercept = kinship_threshold
  ) +
  scale_x_continuous(
    name = "IBS0 [Porportion of SNPs with zero IBS]"
  ) +
  scale_y_continuous(
    name = "Kinship"
  ) +
  theme(
    ggside.panel.scale = 0.15,
    ggside.axis.ticks = element_blank(),
    ggside.axis.text = element_blank(),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank(),
    panel.border = element_blank(),
    ggside.panel.spacing = unit(0, "pt")
  )

plot_name <- "kinship_plot.png"

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

# Kinship plots

write(
  x = "## Relatedness",
  file = md_file,
  append = T
)

write(
  x = "- Pairwise kinship",
  file = md_file,
  append = T
)

density_plot <- ggplot() +
  geom_density(
    data = genomic_relatedness_table,
    mapping = aes(
      y = Kinship
    ),
    fill = "grey90"
  ) +
  geom_vline(
    xintercept = kinship_threshold
  ) +
  scale_x_continuous(
    name = "Kinship"
  ) +
  scale_y_continuous(
    name = "Density"
  ) +
  theme(
    panel.border = element_blank()
  )

plot_name <- "kinship_density.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(density_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

write(
  x = "- Cummulative positive kinship",
  file = md_file,
  append = T
)

cummulative_relatedness_table <- data.frame(
  id = c(genomic_relatedness_table$ID1, genomic_relatedness_table$ID2),
  kinship = c(genomic_relatedness_table$Kinship, genomic_relatedness_table$Kinship)
) %>% 
  mutate(
    kinship = ifelse(kinship > 0, kinship, 0)
  ) %>% 
  summarise(
    cumulated_kinship = sum(kinship),
    .by = "id"
  ) %>% 
  mutate(
    cumulated_kinship = cumulated_kinship / length(ids)
  )

density_plot <- ggplot() +
  geom_density(
    data = cummulative_relatedness_table,
    mapping = aes(
      y = cumulated_kinship
    ),
    fill = "grey90"
  ) +
  geom_vline(
    xintercept = accumulated_kinship_threshold
  ) +
  scale_x_continuous(
    name = "Cumulative kinship"
  ) +
  scale_y_continuous(
    name = "Density"
  ) +
  theme(
    panel.border = element_blank()
  )

plot_name <- "cumulated_kinship_density.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(density_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)


# Remove related individuals

related_ids_table <- genomic_relatedness_table %>% 
  filter(
    Kinship >= kinship_threshold
  ) %>% 
  select(
    ID1, ID2
  ) %>% 
  distinct()

excluded_ids <- c()

relatedness_graph <- graph_from_data_frame(
  related_ids_table,
  directed = F
)

relatedness_components <- components(relatedness_graph)

iteration <- 1

while (length(V(relatedness_graph)) > relatedness_components$no) {
  
  for (component_i in 1:relatedness_components$no) {
    
    print(glue("{Sys.time()}    Percolating relatedness graph - iteration {iteration} ({length(V(relatedness_graph))} nodes in {relatedness_components$no} components)"))
    
    iteration <- iteration + 1
    
    component_nodes <- names(relatedness_components$membership)[relatedness_components$membership == component_i]
    
    if (length(component_nodes) > 1) {
    
    component_table <- related_ids %>% 
      filter(
        ID1 %in% component_nodes & ID2 %in% component_nodes
      )
    
    component_graph <- graph_from_data_frame(
      component_table,
      directed = F
    )
    
    component_components <- components(component_graph)
    
    while (component_components$no == 1) {
    
    component_degree <- degree(component_graph)
    
    most_connected_id <- names(which.max(component_degree))
    
    excluded_ids <- c(excluded_ids, most_connected_id)
    
    component_table <- component_table %>% 
      filter(
        ID1 != most_connected_id & ID2 != most_connected_id
      )
    
    component_graph <- graph_from_data_frame(
      component_table,
      directed = F
    )
    
    component_components <- components(component_graph)
    
    }
    }
  }
  
  related_ids_table <- related_ids_table %>% 
    filter(
      !ID1 %in% excluded_ids & !ID2 %in% excluded_ids
    )
  
  relatedness_graph <- graph_from_data_frame(
    related_ids_table,
    directed = F
  )
  
  relatedness_components <- components(relatedness_graph)
  
}

unrelated_ids <- ids[!ids %in% excluded_ids]

write(
  x = "Percolation of the relatedness graph using a Kinship threshold of {kinship_threshold}: {length(excluded_ids)} excluded, {length(unrelated_ids)} remaining.",
  file = md_file,
  append = T
)


# Compute cumulative kinship among the remaining samples

unrelated_genomic_relatedness_table <- genomic_relatedness_table %>% 
  filter(
    !ID1 %in% excluded_ids & !ID2 %in% excluded_ids
  )

cummulative_relatedness_table <- data.frame(
  id = c(unrelated_genomic_relatedness_table$ID1, unrelated_genomic_relatedness_table$ID2),
  kinship = c(unrelated_genomic_relatedness_table$Kinship, unrelated_genomic_relatedness_table$Kinship)
) %>% 
  mutate(
    kinship = ifelse(kinship > 0, kinship, 0)
  ) %>% 
  summarise(
    cumulated_kinship = sum(kinship),
    .by = "id"
  ) %>% 
  mutate(
    cumulated_kinship = cumulated_kinship / length(ids)
  )

excluded_cumulative_kinship <- cummulative_relatedness_table$id[cummulative_relatedness_table$cumulated_kinship >= accumulated_kinship_threshold]
  

# Kinship plots

print(glue("{Sys.time()}    Relatedness analysis - Writing documentation"))

write(
  x = "## Relatedness after relatedness filtering",
  file = md_file,
  append = T
)

write(
  x = "- Pairwise kinship",
  file = md_file,
  append = T
)

density_plot <- ggplot() +
  geom_density(
    data = unrelated_genomic_relatedness_table,
    mapping = aes(
      y = Kinship
    ),
    fill = "grey90"
  ) +
  geom_vline(
    xintercept = kinship_threshold
  ) +
  scale_x_continuous(
    name = "Kinship"
  ) +
  scale_y_continuous(
    name = "Density"
  ) +
  theme(
    panel.border = element_blank()
  )

plot_name <- "kinship_density_unrelated.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(density_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)

write(
  x = "- Cummulative positive kinship",
  file = md_file,
  append = T
)

density_plot <- ggplot() +
  geom_density(
    data = cummulative_relatedness_table,
    mapping = aes(
      y = cumulated_kinship
    ),
    fill = "grey90"
  ) +
  geom_vline(
    xintercept = accumulated_kinship_threshold
  ) +
  scale_x_continuous(
    name = "Cumulative kinship"
  ) +
  scale_y_continuous(
    name = "Density"
  ) +
  theme(
    panel.border = element_blank()
  )

plot_name <- "cumulated_kinship_density_unrelated.png"

png(
  filename = file.path(docs_dir, plot_name),
  width = 800,
  height = 600
)
grid.draw(density_plot)
device <- dev.off()

write(
  x = paste0("![](", docs_dir_name, "/", plot_name, ")"),
  file = md_file,
  append = T
)


# Get list of samples to exclude

retained_ids <- unrelated_ids[!unrelated_ids %in% excluded_cumulative_kinship]

write(
  x = "Removal of samples with accumulated kinship using threshold of {accumulated_kinship_threshold}: {length(excluded_cumulative_kinship)} excluded, {length(retained_ids)} remaining.",
  file = md_file,
  append = T
)

# Export ids of samples to exclude

to_remove_fam <- fam_data[fam_data[, 2] %in% excluded_ids | fam_data[, 2] %in% excluded_cumulative_kinship, ]

write.table(
  x = to_remove_fam,
  file = exclusion_file,
  append = F,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)

