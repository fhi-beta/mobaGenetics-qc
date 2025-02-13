
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

set.seed(20240112)

debug <- F

if (debug) {
  
  args <- c(
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.12.03/mod8-release_annotation/mod8_pca_both.pcs",
    "/mnt/archive/snpQc/1000Genomes/all_phase3.psam",
    "/mnt/work/oystein/tmp/pca_1kg_moba.md",
    "\"Principal Component Analysys (PCA) vs. 1 KG\"",
    "/mnt/work/oystein/tmp/clusters",
    "/mnt/work/oystein/tmp/ceu_core_ids",
    "/mnt/archive/snpQc/phenotypes/ids_24.08.07.gz",
    "/mnt/archive/moba_genotypes_releases/2024.12.03/batch/moba_genotypes_2024.12.03_batches",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.12.03/mod8-release_annotation/mod8_psam_reconstruction.psam",
    2
  )
  
} else {
  args <- commandArgs(TRUE)
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

plot_folder <- file.path(docs_folder, "plot")

if (!dir.exists(plot_folder)) {
  
  dir.create(plot_folder)
  
}

md_title <- args[4]

cluster_file <- args[5]

ceu_ids_file <- args[6]

id_file <- args[7]

batches_file <- args[8]

psam_file <- args[9]

plink_version <- args[10]


# Libraries
library(janitor)
library(glue)
library(tidyr)
library(dplyr)
library(e1071)
library(ggplot2)
library(ggside)
library(scico)
library(grid)


# Parameters

core_threshold_p <- 0.9


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

# het <- read.table(
#   het_file, 
#   header = F, 
#   col.names = c("iid", "o_hom", "e_hom", "obs_ct", "f")
#   )

id_data  <- read.table(
  file = id_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

batches_data  <- read.table(
  file = batches_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

if(plink_version == 2){
 psam_data  <- read.table(
  file = psam_file,
  header = F,
  sep = "\t",
  col.names = c("fid", "iid", "pat", "mat", "sex"),
  stringsAsFactors = F
)
} else if(plink_version == 1){
  psam_data  <- read.table(
  file = psam_file,
  header = F,
  col.names = c("fid", "iid", "pat", "mat", "sex", "phen"),
  stringsAsFactors = F
)
}





batches_data$batch <- as.factor(batches_data$batch)

populations_order <- c(sort(unique(thousand_genomes_populations$super_pop)), "MoBa")

merged_pcs <- pcs %>% 
  left_join(
    thousand_genomes_populations %>% 
      select(
        iid, pop = super_pop, sub_pop = population
      ),
    by = "iid"
  ) %>% 
  mutate(
    pop = ifelse(is.na(pop), "MoBa", pop),
    pop_factor = factor(pop, levels = populations_order)
  ) %>%
  left_join(
    id_data %>% 
      select(
        iid = sentrix_id, role
      ),
    by = "iid"
  ) %>%
  left_join(
    batches_data,
    by = "iid"
  ) %>%
  left_join(
    psam_data %>% 
      select(
        iid, pat, mat
      ),
    by = "iid"
  ) %>%
  arrange(
    desc(pop_factor)
  )

trios <- subset(merged_pcs, !is.na(pat) & ! is.na(mat))

num_pcs <- 10

pc_child_names <- paste0("pc", 1:num_pcs, "_child")
pc_father_names <- paste0("pc", 1:num_pcs, "_father")
pc_mother_names <- paste0("pc", 1:num_pcs, "_mother")

child_data <- trios %>%
  rename_with(~ pc_child_names, starts_with("pc")) %>%
  select(iid, pat, mat, all_of(pc_child_names), pop_factor_child = pop_factor)

father_data <- subset(merged_pcs, iid %in% trios$pat) %>%
  rename_with(~ pc_father_names, starts_with("pc")) %>%
  select(pat = iid, all_of(pc_father_names), pop_factor_father = pop_factor)


mother_data <- subset(merged_pcs, iid %in% trios$mat) %>%
  rename_with(~ pc_mother_names, starts_with("pc")) %>%
  select(mat = iid, all_of(pc_mother_names), pop_factor_mother = pop_factor)

trios_plot_data <- child_data %>%
  left_join(father_data, by = "pat") %>%
  left_join(mother_data, by = "mat")


write(
  x = paste0("Principal component analysis of the MoBa samples merged with the thousand genomes."),
  file = md_file,
  append = F
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

write(
  x = paste0("| Total | ", nrow(merged_pcs), " |"),
  file = md_file,
  append = T
)


# Plot the PCs against each other

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
  
  moba_data <- subset(merged_pcs, startsWith(pop, "MoBa"))
  
  kg_data <- subset(merged_pcs, !startsWith(pop, "MoBa"))
  
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
  
  print(paste0("Plotting to ", plot_folder, "/", file_name))
  
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

for (i in 1:num_pcs) {
  pc_father <- paste0("pc", i, "_father")
  pc_mother <- paste0("pc", i, "_mother")
  midpoint <- paste0("midpoint_pc", i)
  trios_plot_data[[midpoint]] <- (trios_plot_data[[pc_father]] + trios_plot_data[[pc_mother]]) / 2
}

# Plot the pc of the children against the midpoint between parents, to check that the points cluster along the diagonal:
write(
  x = "## Midpoint of parents vs actual position of children on the PCs (full trios only)",
  file = md_file,
  append = T
)
for (pc_i in 1:num_pcs){
  pc_name <- paste0("pc", pc_i)
  pc_father <- paste0(pc_name, "_father")
  pc_mother <- paste0(pc_name, "_mother")
  midpoint_name <- paste0("midpoint_pc", pc_i)
  child_pc_name <- paste0("pc", pc_i, "_child")
  trios_plot_data[[midpoint_name]] <- (trios_plot_data[[pc_father]] + trios_plot_data[[pc_mother]]) / 2
  trios_plot_data$x <- trios_plot_data[[midpoint_name]]
  trios_plot_data$y <- trios_plot_data[[child_pc_name]]

  plot <- ggplot() +
    theme_bw(
      base_size = 24
    ) +
    
    geom_point(
      data = trios_plot_data,
      mapping = aes(
        x = x,
        y = y
      ),
      alpha = 0.5
    ) +
    
    scale_x_continuous(
      name = paste0("Parents' midpoint on ", pc_name)
    ) +
    scale_y_continuous(
      name = paste0("Child's position on ", pc_name)
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
  
  file_name <- paste0("expected_midpoint_pc_", pc_i, ".png")
  
  print(paste0("Plotting to ", plot_folder, "/", file_name))
  
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



# plot_trios <-function(plot_data, plot_suffix, top_pc){
#  for (pc_i in 1:(top_pc-1)){
#   pc_name_x <- paste0("pc", pc_i)
#   pc_name_y <- paste0("pc", pc_i + 1)

#   child_pc_name_x <- paste0("pc", pc_i, "_child")
#   child_pc_name_y <- paste0("pc", pc_i+1, "_child")

#   father_pc_name_x <- paste0("pc", pc_i, "_father")
#   father_pc_name_y <- paste0("pc", pc_i+1, "_father")

  
#   mother_pc_name_x <- paste0("pc", pc_i, "_mother")
#   mother_pc_name_y <- paste0("pc", pc_i+1, "_mother")

  
#   plot_data$child_x <- plot_data[[child_pc_name_x]]
#   plot_data$child_y <- plot_data[[child_pc_name_y]]
#   plot_data$father_x <- plot_data[[father_pc_name_x]]
#   plot_data$father_y <- plot_data[[father_pc_name_y]]
#   plot_data$mother_x <- plot_data[[mother_pc_name_x]]
#   plot_data$mother_y <- plot_data[[mother_pc_name_y]]



# plot <- ggplot() +
#     theme_bw(
#       base_size = 24
#     ) +
    
#     geom_point(
#       data = plot_data,
#       mapping = aes(
#         x = father_x,
#         y = father_y
#       ),
#       alpha = 0.3,
#       color = "red"
#     ) +
#     geom_point(
#       data = plot_data,
#       mapping = aes(
#         x = mother_x,
#         y = mother_y
#       ),
#       alpha = 0.3,
#       color = "blue"
#     ) +
    
#      geom_segment(data = plot_data, aes(x = child_x, y = child_y, xend = father_x, yend = father_y), linetype = "dashed", color = "red", alpha = 0.2) +
#     geom_segment(data = plot_data, aes(x = child_x, y = child_y, xend = mother_x, yend = mother_y), linetype = "dashed", color = "blue", alpha = 0.2) +
#     geom_point(
#       data = plot_data,
#       mapping = aes(
#         x = child_x,
#         y = child_y
#       ),
#       alpha = 0.5
#     ) +
    
#     scale_x_continuous(
#       name = pc_name_x
#     ) +
#     scale_y_continuous(
#       name = pc_name_y
#     ) +
#     theme(
#       ggside.panel.scale = 0.15,
#       ggside.axis.ticks = element_blank(),
#       ggside.axis.text = element_blank(),
#       ggside.panel.grid = element_blank(),
#       ggside.panel.background = element_blank(),
#       ggside.panel.spacing = unit(0, "pt"),
#       panel.border = element_blank()
#     )
  
#   file_name <- paste0("children_parents_pc", pc_i, "_pc", pc_i + 1, "_", plot_suffix, ".png")
  
#    print(paste0("Plotting to ", plot_folder, "/", file_name))
  
#   png(
#     filename = file.path(plot_folder, file_name),
#     width = 800,
#     height = 600
#   )
#   grid.draw(plot)
#   device <- dev.off()
#   write(
#     x = paste0("![](plot/", file_name, ")"),
#     file = md_file,
#     append = T
#   )

# }
#}


plot_discrete <- function(column, plot_data, top_pc, file_suffix){

 for (pc_i in 1:(top_pc-1)) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  moba_plot_data <- subset(plot_data, startsWith(pop, "MoBa"))
  moba_plot_data$x <- moba_plot_data[[pc_name_x]]
  moba_plot_data$y <- moba_plot_data[[pc_name_y]]
  
  
  kg_plot_data <- subset(merged_pcs, !startsWith(pop, "MoBa"))
  
  kg_plot_data$x <- kg_plot_data[[pc_name_x]]
  kg_plot_data$y <- kg_plot_data[[pc_name_y]]
  
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
      data = moba_plot_data,
      mapping = aes(
        x = x,
        y = y,
        col = .data[[column]]
      ),
      alpha = 0.2
    ) +
    geom_density2d(
      data = kg_plot_data,
      mapping = aes(
       x = x,
       y = y
      ),
      alpha = 0.8
    ) +
     geom_xsidedensity(
       data = moba_plot_data,
       mapping = aes(
         x = x,
         y = after_stat(density),
         fill = .data[[column]]
       ),
       alpha = 0.8
     ) +
     geom_ysidedensity(
       data = moba_plot_data,
       mapping = aes(
         x = after_stat(density),
         y = y,
         fill = .data[[column]]
       ),
       alpha = 0.8
     ) +
    scale_x_continuous(,
      name = pc_name_x
    ) +
    scale_y_continuous(
      name = pc_name_y
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
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_", file_suffix, ".png")
  
  print(paste0("Plotting to ", plot_folder, "/", file_name))
  
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

}

#plot_trios(trios_plot_data, "trios", 10)

#plot_discrete("pop_factor", merged_pcs, 3, "pop_factor")

write(
  x = "## PCs with batches marked",
  file = md_file,
  append = T
)

plot_discrete("batch", merged_pcs, 3, "batch")

# 1kg cluster size
kg <- subset(merged_pcs, !startsWith(pop, "MoBa"))

kg <- kg %>%  
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

train_df <- as.data.frame(kg_1) %>% 
  mutate(
    pop_factor = factor(pop_kg)
  ) %>% 
  select(
    pop_factor, starts_with("pc")
  )

classifier <- svm(
  formula = pop_factor ~ ., 
  data = train_df, 
  type = 'C-classification', 
  kernel = 'linear', 
  probability = TRUE
)

test_df <- as.data.frame(kg_2) %>% 
  select(
    starts_with("pc")
  )

pop_inference <- predict(
  classifier, 
  newdata = test_df, 
  probability = TRUE
  )

kg_2$pop_inference <- as.character(pop_inference)

probabilities <- as.data.frame(attr(pop_inference, "probabilities"))
names(probabilities) <- paste0("p_", names(probabilities))

kg_2 <- cbind(kg_2, probabilities)

write(
  x = paste0("### Clustering in the 1KG"),
  file = md_file,
  append = T
)

write(
  x = paste0("Clustering of the 1KG after an 80-20 split using svm in the top 10 PCs."),
  file = md_file,
  append = T
)

pop_table <- as.data.frame(table(kg_2$pop_kg, kg_2$pop_inference))
pop_table$Var1 <- as.character(pop_table$Var1)
pop_table$Var2 <- as.character(pop_table$Var2)
populations <- sort(unique(c(pop_table$Var1, pop_table$Var2)))

write(
  x = paste0("| Population | ", paste(populations, collapse = " | "), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| - |", paste(rep(" - ", length(populations)), collapse = " | "), " |"),
  file = md_file,
  append = T
)

for (population_ref in populations) {
  
  line <- paste0("| ", population_ref)
  
  for (population_inferred in populations) {
    
    n <- sum(kg_2$pop_kg == population_ref & kg_2$pop_inference == population_inferred)
    
    line <- paste0(line, " | ", n)
    
  }
  
  line <- paste0(line, " |")
  
  write(
    x = line,
    file = md_file,
    append = T
  )
  
}

write(
  x = "\n",
  file = md_file,
  append = T
)


kg_population_probabilities_plot <- kg_2 %>% 
  select(
    pop_kg, starts_with("p_")
    ) %>% 
  pivot_longer(
    cols = starts_with("p_"),
    names_prefix = "p_",
    names_to = "pop_pred",
    values_to = "probability"
  ) %>% 
  mutate(
    pop_kg_factor = factor(pop_kg),
    pop_pred_factor = factor(pop_pred)
  )

kg_pop_plot <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_violin(
    data = kg_population_probabilities_plot,
    mapping = aes(
      x = probability,
      y = pop_pred_factor,
      col = pop_pred_factor
    )
  ) +
  scale_x_continuous(
    name = "Probability"
  ) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) +
  facet_grid(
    pop_kg_factor ~ .
  )

file_name <- "kg_pop_plot.png"

png(
  filename = file.path(plot_folder, file_name),
  width = 800,
  height = 600
)
grid.draw(kg_pop_plot)
device <- dev.off()

write(
  x = paste0("![](plot/", file_name, ")"),
  file = md_file,
  append = T
)


# Inference in MoBa
train_df <- subset(merged_pcs, !startsWith(pop, "MoBa"))

train_df <- train_df %>%  
  mutate(
    pop_factor = factor(pop)
  ) %>% 
  select(
    starts_with("pc"),
    pop_factor
  )

classifier <- svm(
  formula = pop_factor ~ ., 
  data = train_df, 
  type = 'C-classification', 
  kernel = 'linear', 
  probability = TRUE
)

moba_df <- subset(merged_pcs, startsWith(pop, "MoBa")) %>% 
  select(
    fid, iid, starts_with("pc"), role
  )

predict_df <- moba_df %>% 
  select(
    starts_with("pc")
  )

pop_inference <- predict(
  classifier, 
  newdata = predict_df, 
  probability = TRUE
)

moba_df$pop_inference <- as.character(pop_inference)

probabilities <- as.data.frame(attr(pop_inference, "probabilities"))
names(probabilities) <- paste0("p_", names(probabilities))

moba_df <- cbind(moba_df, probabilities)

moba_df$pop_inference[moba_df$pop_inference == "EUR" & moba_df$p_EUR >= core_threshold_p] <- "EUR_core"


write(
  x = paste0("### Clustering in MoBa"),
  file = md_file,
  append = T
)

write(
  x = paste0("Clustering of MoBa participants with svm using the top 10 PCs and 1KG as training."),
  file = md_file,
  append = T
)

pop_table <- as.data.frame(table(moba_df$pop_inference))
pop_table$Var1 <- as.character(pop_table$Var1)
populations <- sort(pop_table$Var1)

write(
  x = paste0("| Population | ", paste(populations, collapse = " | "), " |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| - |", paste(rep(" - ", length(populations)), collapse = " | "), " |"),
  file = md_file,
  append = T
)
  
  line <- paste0("| MoBa")
  
  for (population_inferred in populations) {
    
    n <- sum(moba_df$pop_inference == population_inferred)
    
    line <- paste0(line, " | ", n)
    
  }
  
  line <- paste0(line, " |")
  
  write(
    x = line,
    file = md_file,
    append = T
  )
  
  write(
    x = "\n",
    file = md_file,
    append = T
  )


moba_population_probabilities_plot <- moba_df %>% 
  select(
    pop_inference, starts_with("p_")
  ) %>% 
  pivot_longer(
    cols = starts_with("p_"),
    names_prefix = "p_",
    names_to = "pop_pred",
    values_to = "probability"
  ) %>% 
  mutate(
    pop_inference_factor = factor(pop_inference),
    pop_pred_factor = factor(pop_pred)
  )

moba_pop_plot <- ggplot() +
  theme_bw(
    base_size = 24
  ) +
  geom_violin(
    data = moba_population_probabilities_plot,
    mapping = aes(
      x = probability,
      y = pop_pred_factor,
      col = pop_pred_factor
    )
  ) +
  scale_x_continuous(
    name = "Probability"
  ) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) +
  facet_grid(
    pop_inference_factor ~ .
  )

file_name <- "moba_pop_plot.png"

png(
  filename = file.path(plot_folder, file_name),
  width = 800,
  height = 600
)
grid.draw(moba_pop_plot)
device <- dev.off()

write(
  x = paste0("![](plot/", file_name, ")"),
  file = md_file,
  append = T
)


populations_colors <- scico(
  n = length(populations_order),
  begin = 0.2,
  end = 0.8,
  palette = "hawaii"
)

write(
  x = paste0("### Principal components with F statistics"),
  file = md_file,
  append = T
)

plot_folder <- file.path(docs_folder, "plot")










# Plot the PCs of the clusters



write(
  x = paste0("### Principal components with color of the top cluster"),
  file = md_file,
  append = T
)

plot_folder <- file.path(docs_folder, "plot")

populations <- sort(unique(c(moba_df$pop_inference, merged_pcs$pop)))

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  moba_plot_data <- moba_df
  moba_plot_data$x <- moba_df[[pc_name_x]]
  moba_plot_data$y <- moba_df[[pc_name_y]]
  moba_plot_data$pop_factor <- factor(moba_plot_data$pop, levels = populations)
  
  kg_plot_data <- subset(merged_pcs, !startsWith(pop, "MoBa")) %>% 
    mutate(
      pop_factor = factor(pop, levels = populations)
    )
  
  kg_plot_data$x <- kg_plot_data[[pc_name_x]]
  kg_plot_data$y <- kg_plot_data[[pc_name_y]]
  
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
      data = moba_plot_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop_factor
      ),
      alpha = 0.1
    ) +
    geom_density2d(
      data = kg_plot_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop_factor
      )
    ) +
    geom_xsidedensity(
      data = kg_plot_data,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = kg_plot_data,
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
  
  print(paste0("Plotting to ", plot_folder, "/", file_name))
  
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
  moba_df,
  file = destination_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)

core_ids <- moba_df %>% 
  filter(
    pop_inference == "EUR_core"
  ) %>% 
  select(
    fid, iid
  )


very_mysterious_samples <- merged_pcs %>% 
  filter(
    pop == "MoBa_Very_Mysterious"
  ) %>% 
  select(
    fid, iid
  )

write.table(
  core_ids,
  file = ceu_ids_file,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)

write.table(
  very_mysterious_samples,
  file = "/mnt/work/oystein/tmp/very_mysterious_samples.txt",
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)