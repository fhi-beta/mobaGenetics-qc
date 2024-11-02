
##
#
# This script plots the PCA of MoBa with batches annotated in color.
#
##

# Parameters
md_title <- "Principal Component Analysys (PCA) vs. batch"

# Files
moba_pcs_file <- "/mnt/archive/moba_genotypes_releases/2024.07.01/pca/moba_genotypes_2024.07.01.pcs"
moba_batches_file <- "/mnt/archive/moba_genotypes_releases/2024.07.01/batch/moba_genotypes_2024.07.01_batches"
md_file <- "/mnt/work/marc/github/tmp/mobaGenetics-qc/qc-pipeline/docs/2024.07.01/pca_batches.md"

docs_folder <- dirname(md_file)

if (!dir.exists(docs_folder)) {

  dir.create(docs_folder)

}

plot_folder <- file.path(docs_folder, "plot")

if (!dir.exists(plot_folder)) {

  dir.create(plot_folder)

}

# Libraries
library(janitor)
library(glue)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggside)
library(scico)
library(grid)


# load data
moba_pcs <- read.table(
  file = moba_pcs_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>%
  clean_names() %>%
  filter(
    fid != 0
  )

moba_pcs <- moba_pcs[sample(1:nrow(moba_pcs), 10000), ]

moba_batches <- read.table(
  file = moba_batches_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>%
  clean_names()

merged_pcs <- moba_pcs %>%
    left_join(
        moba_batches,
        by = "iid"
    )

batches_order <- unique(merged_pcs$batch)

# Write docs

write(
  x = paste0("# ", md_title),
  file = md_file,
  append = F
)

write(
  x = paste0("Principal component analysis of the MoBa samples vs genotyping batch."),
  file = md_file,
  append = T
)

write(
  x = paste0("| Batch | Number of samples |"),
  file = md_file,
  append = T
)
write(
  x = paste0("| ----- | ----------------- |"),
  file = md_file,
  append = T
)

n_samples <- as.data.frame(table(merged_pcs$batch))

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

batches_colors <- scico(
  n = length(batches_order),
  begin = 0.2,
  end = 0.8,
  palette = "batlow"
)


for (pc_i in 1:9) {

  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)

  merged_pcs$x <- merged_pcs[[pc_name_x]]
  merged_pcs$y <- merged_pcs[[pc_name_y]]

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
      data = merged_pcs,
      mapping = aes(
        x = x,
        y = y,
        col = batch
      ),
      alpha = 0.1
    ) +
    geom_xsidedensity(
      data = merged_pcs,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = batch
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = merged_pcs,
      mapping = aes(
        x = after_stat(density),
        y = y,
        fill = batch
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
      values = batches_colors,
      drop = F
    ) +
    scale_fill_manual(
      name = "Population",
      values = batches_colors,
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

  file_name <- paste0(pc_name_x, "_", pc_name_y, "_batch.png")

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




