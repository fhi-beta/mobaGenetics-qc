
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

args <- commandArgs(TRUE)

if (length(args) != 4) {
  
  stop(paste0("Four arguments expected: pcs file, thousand genomes population file, md file, md title. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

pcs_file <- args[1]

if (!file.exists(pcs_file)) {
  
  stop("PCs file not found")
  
}

thousand_genomes_populations_file <- args[2]

if (!file.exists(thousand_genomes_populations_file)) {
  
  stop("Thousand genomes population file not found")
  
}

md_file <- args[3]
docs_folder <- dirname(md_file)

if (!dir.exists(docs_folder)) {
  
  dir.create(docs_folder)
  
}

md_title <- args[4]

cluster_file <- args[5]

ceu_ids_file <- args[6]

n_pcs <- args[7]

if (n_pcs < 1 | n_pcs > 10) {
  
  stop(paste0("The number of PCs should be between 1 and 10. Input: ", n_pcs, "."))
  
}


# Local debug - do not uncomment
# 
# pcs_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-population-clustering/snp012/pca_both.pcs"
# thousand_genomes_populations_file <- "/mnt/archive/snpQc/1000Genomes/all_phase3.psam"
# md_file <- "/mnt/work/marc/github/mobaGenetics-qc/qc-pipeline/docs/snp012/pca_1kg_moba.md"
# md_title <- "Principal Component Analysys (PCA) in snp012 vs. 1kg"


# Libraries
library(janitor)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggside)
library(scico)
library(grid)


# Load data
thousand_genomes_populations <- read.table(
  file = thousand_genomes_populations_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  rename(
    iid = x_iid
  )

pcs <- read.table(
  file = pcs_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names()


# Merge

populations_order <- c(sort(unique(thousand_genomes_populations$super_pop)), "MoBa")

merged_pcs <- pcs %>% 
  left_join(
    thousand_genomes_populations %>% 
      select(
        iid, pop = super_pop
      ),
    by = "iid"
  ) %>% 
  mutate(
    pop = ifelse(is.na(pop), "MoBa", pop),
    pop_factor = factor(pop, levels = populations_order)
  ) %>% 
  arrange(
    desc(pop_factor)
  )


# Write docs

write(
  x = paste0("# ", md_title),
  file = md_file,
  append = F
)

write(
  x = paste0("Principal component analysis of the MoBa samples merged with the thousand genomes."),
  file = md_file,
  append = T
)

write(
  x = paste0("| Population code | Description |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| --------------- | ----------- |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| AFR | African |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| AMR | Admixed American |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| EAS | East Asian |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| EUR | European |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| MID | Middle Eastern |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| SAS | South Asian |"),
  file = md_file,
  append = T
)
write(
  x = paste0("\n\n"),
  file = md_file,
  append = T
)



write(
  x = paste0("| Population | Number of samples |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| ---------- | ----------------- |"),
  file = md_file,
  append = T
)

n_samples <- as.data.frame(table(merged_pcs$pop_factor))

for (sample_i in 1:nrow(n_samples)) {
  
  write(
    x = paste0("| ", n_samples[sample_i, 1], " | ", n_samples[sample_i, 2], " |"),
    file = md_file,
    append = T
  )
  
}


# Plot the PCs against each other

plot_folder <- file.path(docs_folder, "plot")

if (!dir.exists(plot_folder)) {
  
  dir.create(plot_folder)
  
}

kg_populations_colors <- scico(
  n = length(populations_order) - 1,
  begin = 0.2,
  end = 0.8,
  palette = "hawaii"
)

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  merged_pcs$x <- merged_pcs[[pc_name_x]]
  merged_pcs$y <- merged_pcs[[pc_name_y]]
  
  moba_data <- merged_pcs %>% 
    filter(
      pop == "MoBa"
    )
  
  kg_data <- merged_pcs %>% 
    filter(
      pop != "MoBa"
    )
  
  write(
    x = paste0("### ", pc_name_y, " vs. ", pc_name_x),
    file = md_file,
    append = T
  )
  
  plot <- ggplot() +
    theme_bw(
      base_size = 24
    ) +
    geom_point(
      data = moba_data,
      mapping = aes(
        x = x,
        y = y
      ),
      alpha = 0.1,
        col = "black"
    ) +
    geom_density2d(
      data = kg_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop_factor
      )
    ) +
    geom_xsidedensity(
      data = merged_pcs,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = merged_pcs,
      mapping = aes(
        x = after_stat(density),
        y = y,
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    scale_x_continuous(
      name = pc_name_x
    ) +
    scale_y_continuous(
      name = pc_name_y
    ) +
    scale_color_manual(
      name = "Population",
      values = c(kg_populations_colors, "black"),
      drop = F
    ) +
    scale_fill_manual(
      name = "Population",
      values = c(kg_populations_colors, "grey80"),
      drop = F
    ) +
    theme(
      ggside.panel.scale = 0.15,
      ggside.axis.ticks = element_blank(),
      ggside.axis.text = element_blank(),
      ggside.panel.grid = element_blank(),
      ggside.panel.background = element_blank(),
      ggside.panel.spacing = unit(0, "pt"),
      panel.border = element_blank()
    )
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_1kg.png")
  
  print(paste0("Plotting to ", plot_folder, file_name))
  
  png(
    filename = file.path(plot_folder, file_name),
    width = 800,
    height = 600
  )
  grid.draw(plot)
  device <- dev.off()
  
  write(
    x = paste0("![](plot/", file_name, ")"),
    file = md_file,
    append = T
  )
  
}


# 1kg cluster size

kg <- merged_pcs %>% 
  filter(
    pop != "MoBa"
  ) %>% 
  select(
    iid_kg = iid,
    pc1_kg = pc1,
    pc2_kg = pc2,
    pc3_kg = pc3,
    pc4_kg = pc4,
    pc5_kg = pc5,
    pc6_kg = pc6,
    pc7_kg = pc7,
    pc8_kg = pc8,
    pc9_kg = pc9,
    pc10_kg = pc10,
    pop_kg = pop
  ) %>% 
  mutate(
    merge = 0
  )

kg_1 <- kg %>% 
  group_by(
    pop_kg
  ) %>% 
  sample_frac(
    0.2
  )

kg_2 <- kg %>% 
  filter(
    !iid_kg %in% kg_1$iid_kg
  )

names(kg_1) <- paste0(names(kg_1), "_1")
names(kg_2) <- paste0(names(kg_2), "_2")

names(kg_1)[ncol(kg_1)] <- "merge"
names(kg_2)[ncol(kg_2)] <- "merge"

kg_distance <- kg_1 %>% 
  full_join(
    kg_2,
    by = "merge",
    multiple = "all"
  ) %>% 
  select(
    -merge
  ) %>% 
  mutate(
    distance = sqrt(
      (pc1_kg_1 - pc1_kg_2)^2 +
        (pc2_kg_1 - pc2_kg_2)^2 +
        (pc3_kg_1 - pc3_kg_2)^2 +
        (pc4_kg_1 - pc4_kg_2)^2 +
        (pc5_kg_1 - pc5_kg_2)^2 +
        (pc6_kg_1 - pc6_kg_2)^2 +
        (pc7_kg_1 - pc7_kg_2)^2 +
        (pc8_kg_1 - pc8_kg_2)^2 +
        (pc9_kg_1 - pc9_kg_2)^2 +
        (pc10_kg_1 - pc10_kg_2)^2
    )
  )

kg_distance$distance <- 0

for (pc in 1:n_pcs) {
  
  pc_kg_1 <- paste0("pc", pc, "_kg_1")
  pc_kg_2 <- paste0("pc", pc, "_kg_2")
  
  kg_distance$distance <- distance_matrix$distance + ((distance_matrix[[pc_moba]] - distance_matrix[[pc_kg]]) ^ 2)
  
}

kg_distance$distance <- sqrt(kg_distance$distance)

kg_population_distance <- kg_distance %>%
  select(
    pop_kg_1, pop_kg_2, distance
  ) %>% 
  group_by(
    pop_kg_1, pop_kg_2
  ) %>% 
  summarize(
    mean = mean(distance),
    sd = sd(distance),
    q5 = quantile(distance, 0.05),
    q50 = quantile(distance, 0.5),
    q95 = quantile(distance, 0.95),
    .groups = "keep"
  )

write(
  x = paste0("### Cluster distance in the 1KG"),
  file = md_file,
  append = T
)

kg_population_distance_plot <- kg_population_distance %>% 
  mutate(
    pop_kg_1_factor = factor(pop_kg_1),
    pop_kg_2_factor = factor(pop_kg_2)
  )

distance_plot <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_point(
    data = kg_population_distance_plot,
    mapping = aes(
      x = q50,
      y = pop_kg_2_factor,
      col = pop_kg_2_factor
    )
  ) +
  geom_segment(
    data = kg_population_distance_plot,
    mapping = aes(
      x = q5,
      xend = q95,
      y = pop_kg_2_factor,
      yend = pop_kg_2_factor,
      col = pop_kg_2_factor
    )
  ) +
  scale_x_continuous(
    name = "Distance [median - 90%]"
  ) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) +
  facet_grid(
    pop_kg_1_factor ~ .
  )

file_name <- "kg_distance.png"

png(
  filename = file.path(plot_folder, file_name),
  width = 800,
  height = 600
)
grid.draw(distance_plot)
device <- dev.off()

write(
  x = paste0("![](plot/", file_name, ")"),
  file = md_file,
  append = T
)


# Distance between MoBa and 1kg

moba <- merged_pcs %>% 
  filter(
    pop == "MoBa"
  ) %>% 
  select(
    iid_moba = iid,
    pc1_moba = pc1,
    pc2_moba = pc2,
    pc3_moba = pc3,
    pc4_moba = pc4,
    pc5_moba = pc5,
    pc6_moba = pc6,
    pc7_moba = pc7,
    pc8_moba = pc8,
    pc9_moba = pc9,
    pc10_moba = pc10
  ) %>% 
  mutate(
    merge = 0
  )

kg <- merged_pcs %>% 
  filter(
    pop != "MoBa"
  ) %>% 
  select(
    iid_kg = iid,
    pc1_kg = pc1,
    pc2_kg = pc2,
    pc3_kg = pc3,
    pc4_kg = pc4,
    pc5_kg = pc5,
    pc6_kg = pc6,
    pc7_kg = pc7,
    pc8_kg = pc8,
    pc9_kg = pc9,
    pc10_kg = pc10,
    pop_kg = pop
  ) %>% 
  mutate(
    merge = 0
  )

distance_matrix <- moba %>% 
  full_join(
    kg,
    by = "merge",
    multiple = "all"
  ) %>% 
  select(
    -merge
  )

distance_matrix$distance <- 0

for (pc in 1:n_pcs) {
  
  pc_moba <- paste0("pc", pc, "_moba")
  pc_kg <- paste0("pc", pc, "_kg")
  
  distance_matrix$distance <- distance_matrix$distance + ((distance_matrix[[pc_moba]] - distance_matrix[[pc_kg]]) ^ 2)
  
}

distance_matrix$distance <- sqrt(distance_matrix$distance)

population_distance_matrix <- distance_matrix %>% 
  select(
    iid_moba, pop_kg, distance
  ) %>% 
  group_by(
    iid_moba, pop_kg
  ) %>% 
  summarize(
    mean = mean(distance),
    sd = sd(distance),
    q5 = quantile(distance, 0.05),
    q50 = quantile(distance, 0.5),
    q95 = quantile(distance, 0.95),
    .groups = "keep"
  )
  
population_distance_matrix <- population_distance_matrix %>% 
  left_join(
    kg_population_distance %>% 
      ungroup() %>% 
      filter(
        pop_kg_1 == pop_kg_2
      ) %>% 
      select(
        pop_kg = pop_kg_1,
        q95_kg = q95
      ),
    by = "pop_kg"
  )

population_inference <- population_distance_matrix %>% 
  mutate(
    population_cluster = ifelse(q95 <= q95_kg, pop_kg, NA)
  ) %>% 
  filter(
    !is.na(population_cluster)
  ) %>% 
  select(
    iid_moba, population_cluster
  ) %>% 
  arrange(
    population_cluster
    ) %>% 
  group_by(
    iid_moba
  ) %>% 
  summarize(
    population_cluster = paste(population_cluster, collapse = "_"),
    .groups = "keep"
  )

admixed <- data.frame(
  iid_moba = unique(population_distance_matrix$iid_moba[!population_distance_matrix$iid_moba %in% population_inference$iid_moba]),
  population_cluster = "NONE"
)

population_inference <- rbind(population_inference, admixed) %>% 
  rename(
    iid = iid_moba
  )


# Plot the PCs of the clusters

write(
  x = paste0("### Population clustering in MoBa"),
  file = md_file,
  append = T
)

write(
  x = paste0("Clustering using nearest neighbors in the top ", n_pcs, " PCs."),
  file = md_file,
  append = T
)

plot_folder <- file.path(docs_folder, "plot")

kg_populations_colors <- scico(
  n = length(populations_order) - 1,
  begin = 0.2,
  end = 0.8,
  palette = "hawaii"
)

populations_order <- c(populations_order, "NONE")

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  merged_pcs$x <- merged_pcs[[pc_name_x]]
  merged_pcs$y <- merged_pcs[[pc_name_y]]
  
  moba_data <- merged_pcs %>% 
    filter(
      pop == "MoBa"
    ) %>% 
    left_join(
      population_inference,
      by = "iid"
    ) %>% 
    mutate(
      pop_factor = factor(population_cluster, levels = populations_order)
    )
  
  kg_data <- merged_pcs %>% 
    filter(
      pop != "MoBa"
    ) %>% 
    mutate(
      pop_factor = factor(pop, levels = populations_order)
    )
  
  write(
    x = paste0("### ", pc_name_y, " vs. ", pc_name_x),
    file = md_file,
    append = T
  )
  
  plot <- ggplot() +
    theme_bw(
      base_size = 24
    ) +
    geom_point(
      data = moba_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop_factor
      ),
      alpha = 0.1
    ) +
    geom_density2d(
      data = kg_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop_factor
      )
    ) +
    geom_xsidedensity(
      data = merged_pcs,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = merged_pcs,
      mapping = aes(
        x = after_stat(density),
        y = y,
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    scale_x_continuous(
      name = pc_name_x
    ) +
    scale_y_continuous(
      name = pc_name_y
    ) +
    scale_color_manual(
      name = "Population",
      values = c(kg_populations_colors, "black", "grey80"),
      drop = F
    ) +
    scale_fill_manual(
      name = "Population",
      values = c(kg_populations_colors, "grey80", "grey80"),
      drop = F
    ) +
    theme(
      ggside.panel.scale = 0.15,
      ggside.axis.ticks = element_blank(),
      ggside.axis.text = element_blank(),
      ggside.panel.grid = element_blank(),
      ggside.panel.background = element_blank(),
      ggside.panel.spacing = unit(0, "pt"),
      panel.border = element_blank()
    )
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_1kg_inferred.png")
  
  print(paste0("Plotting to ", plot_folder, file_name))
  
  png(
    filename = file.path(plot_folder, file_name),
    width = 800,
    height = 600
  )
  grid.draw(plot)
  device <- dev.off()
  
  write(
    x = paste0("![](plot/", file_name, ")"),
    file = md_file,
    append = T
  )
  
}

  
if (endsWith(cluster_file, ".gz")) {
  
  destination_file <- gzfile(cluster_file)
  
} else {
  
  destination_file <- cluster_file
  
}

write.table(
  population_inference %>% 
    select(
      iid, population_cluster
    ),
  file = destination_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

ceu_ids <- population_inference$iid[population_inference$population_cluster == "CEU"]

ceu_ids_plink <- pcs %>% 
  select(
    fid, iid
  ) %>% 
  filter(
    iid %in% ceu_ids
  )

write.table(
  ceu_ids_plink,
  file = ceu_ids_file,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)
