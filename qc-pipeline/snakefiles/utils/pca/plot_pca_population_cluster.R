
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

args <- commandArgs(TRUE)

if (length(args) != 5) {
  
  stop(paste0("Four arguments expected: pcs file, thousand genomes population file, md file, md title, output file. ", length(args), " found: ", paste(args, collapse = ", ")))
  
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
export_file <- args[5]


# Local debug - do not uncomment
# 
# pcs_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/pca_both.pcs"
# thousand_genomes_populations_file <- "/mnt/archive/snpQc/1000Genomes/all_phase3.psam"
# md_file <- "/mnt/work/marc/github/mobaGenetics-qc/qc-pipeline/docs/snp014/pca_1kg_moba.md"
# md_title <- "Principal Component Analysys (PCA) in snp014"
# export_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/moba.populations"


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

plot_data <- pcs %>% 
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
  x = paste0("| Population | Number of samples |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| ---------- | ----------------- |"),
  file = md_file,
  append = T
)

n_samples <- as.data.frame(table(plot_data$pop_factor))

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
  
  plot_data$x <- plot_data[[pc_name_x]]
  plot_data$y <- plot_data[[pc_name_y]]
  
  moba_data <- plot_data %>% 
    filter(
      pop == "MoBa"
    )
  
  kg_data <- plot_data %>% 
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
      data = plot_data,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = pop_factor
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = plot_data,
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
      values = kg_populations_colors
    ) +
    scale_fill_manual(
      name = "Population",
      values = c(kg_populations_colors, "grey80")
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
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, ".png")
  
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


# Export MoBa populations

moba_populations <- plot_data %>% 
  filter(
    pop == "MoBa"
  ) %>% 
  select(fid, iid, pop)

write.table(
  x = moba_populations,
  file = export_file,
  col.names = T,
  row.names = F,
  sep = "\t"
)
