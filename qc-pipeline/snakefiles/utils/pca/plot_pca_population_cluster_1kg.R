
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

set.seed(20240112)

debug <- T

if (debug) {
  
  args <- c(
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.12.03/mod8-release_annotation/mod8_pca_both.pcs",
    "/mnt/archive/snpQc/1000Genomes/all_phase3.psam",
    "/mnt/work/qc_genotypes/pipeOut_dev/2024.12.03/mod8-release_annotation/mod8_best_snps.het",
    "/mnt/work/oystein/tmp/pca_1kg_moba.md",
    "\"Principal Component Analysys (PCA) vs. 1 KG\"",
    "/mnt/work/oystein/tmp/clusters",
    "/mnt/work/oystein/tmp/ceu_core_ids",
    4
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

het_file <- args[3]

if (!file.exists(het_file)) {
  
  stop("Heterozygosity file not found")
  
}

md_file <- args[4]
docs_folder <- dirname(md_file)

if (!dir.exists(docs_folder)) {
  
  dir.create(docs_folder)
  
}

plot_folder <- file.path(docs_folder, "plot")

if (!dir.exists(plot_folder)) {
  
  dir.create(plot_folder)
  
}

md_title <- args[5]

cluster_file <- args[6]

ceu_ids_file <- args[7]

std_cutoff <- as.numeric(args[8])


# Local debug - do not uncomment
# 
# pcs_file <- "/mnt/archive/snpQc/pipeOut_dev_2024.01.05/mod3-population-clustering/snp012/pca_both.pcs"
# thousand_genomes_populations_file <- "/mnt/archive/snpQc/1000Genomes/all_phase3.psam"
# md_file <- "/mnt/work/marc/github/mobaGenetics-qc/qc-pipeline/docs/snp012/pca_1kg_moba.md"
# md_title <- "Principal Component Analysys (PCA) in snp012 vs. 1kg"


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

het <- read.table(
  het_file, 
  header = F, 
  col.names = c("iid", "o_hom", "e_hom", "obs_ct", "f")
  )

het$het_rate <- (het$obs_ct - het$o_hom)/het$obs_ct
het$stds_het_rate <- 0
mean_het_rate <- mean(het$het_rate)
std_het_rate <- sd(het$het_rate)

for (std in 1:(2*(std_cutoff-1))){
  het$stds_het_rate[het$het_rate > mean_het_rate + 0.5*std*std_het_rate] <- 0.5*std 
}


# Merge

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
    het %>% 
      select(
        iid, het_rate, stds_het_rate, f
      ),
    by = "iid"
  ) %>% 
  arrange(
    desc(pop_factor)
  )
# Testing filtering for < 0.5 std het rate
#merged_pcs <- subset(merged_pcs, (pop_factor == "MoBa" & stds_het_rate == 0) | pop_factor != "MoBa")
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

train_df <- merged_pcs %>% 
  filter(
    pop != "MoBa"
  ) %>% 
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

moba_df <- merged_pcs %>% 
  filter(
    pop == "MoBa"
  ) %>% 
  select(
    fid, iid, starts_with("pc"), stds_het_rate, het_rate, f
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

#testing filtering for f<0 and not EUR_core

#moba_df <- subset(moba_df, f>0 | (f<0 & pop_inference == "EUR_core"))

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


# Plot the F stats rates vs PCs
# TODO: Move this into a function for plotting other attributes against PCs

write(
  x = paste0("### Principal components with F statistics"),
  file = md_file,
  append = T
)

plot_folder <- file.path(docs_folder, "plot")

# stds <- sort(unique(moba_df$stds_het_rate))
# moba_df$stds_het_rate <- factor(moba_df$stds_het_rate, levels = stds)
# populations <- sort(unique(c(moba_df$pop_inference, merged_pcs$pop)))

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  moba_plot_data <- moba_df
  moba_plot_data$x <- moba_df[[pc_name_x]]
  moba_plot_data$y <- moba_df[[pc_name_y]]
  
  
  kg_plot_data <- merged_pcs %>% 
    filter(
      pop != "MoBa"
    ) %>% 
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
        col = f
      ),
      alpha = 0.4
    ) +
    geom_point(
      data = kg_plot_data,
      mapping = aes(
        x = x,
        y = y
      ),
      alpha = 0.5
    ) +
    # geom_xsidedensity(
    #   data = moba_plot_data,
    #   mapping = aes(
    #     x = x,
    #     y = after_stat(density),
    #     fill = stds_het_rate
    #   ),
    #   alpha = 0.8
    # ) +
    # geom_ysidedensity(
    #   data = moba_plot_data,
    #   mapping = aes(
    #     x = after_stat(density),
    #     y = y,
    #     fill = stds_het_rate
    #   ),
    #   alpha = 0.8
    # ) +
    scale_color_stepsn(
      breaks = c(-0.02, 0),
      colours = c("red", "blue", "green")
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
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_1kg_f_stat.png")
  
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


# Plot the heterozygote rates vs PCs
# TODO: Move this into a function for plotting other attributes against PCs

write(
  x = paste0("### Principal components with heterozygote rate"),
  file = md_file,
  append = T
)

plot_folder <- file.path(docs_folder, "plot")

stds <- sort(unique(moba_df$stds_het_rate))
moba_df$stds_het_rate <- factor(moba_df$stds_het_rate, levels = stds)
populations <- sort(unique(c(moba_df$pop_inference, merged_pcs$pop)))

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  moba_plot_data <- moba_df
  moba_plot_data$x <- moba_df[[pc_name_x]]
  moba_plot_data$y <- moba_df[[pc_name_y]]
  
  
  kg_plot_data <- merged_pcs %>% 
    filter(
      pop != "MoBa"
    ) %>% 
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
        col = stds_het_rate
      ),
      alpha = 0.4
    ) +
    geom_point(
      data = kg_plot_data,
      mapping = aes(
        x = x,
        y = y
      ),
      alpha = 0.2
    ) +
    geom_xsidedensity(
      data = moba_plot_data,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = stds_het_rate
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = moba_plot_data,
      mapping = aes(
        x = after_stat(density),
        y = y,
        fill = stds_het_rate
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
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_1kg_heterozygote_stds.png")
  
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
  
  kg_plot_data <- merged_pcs %>% 
    filter(
      pop != "MoBa"
    ) %>% 
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

write.table(
  core_ids,
  file = ceu_ids_file,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)
