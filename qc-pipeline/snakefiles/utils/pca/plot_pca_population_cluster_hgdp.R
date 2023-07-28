
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

args <- commandArgs(TRUE)

# if (length(args) != 5) {
#   
#   stop(paste0("Four arguments expected: pcs file, thousand genomes population file, md file, md title, output file. ", length(args), " found: ", paste(args, collapse = ", ")))
#   
# }

moba_scores_file <- args[1]

if (!file.exists(moba_scores_file)) {
  
  stop("MoBa scores file not found")
  
}

hgdp_scores_file <- args[2]

if (!file.exists(hgdp_scores_file)) {
  
  stop("HGDP scores file not found")
  
}

md_file <- args[3]
docs_folder <- dirname(md_file)

if (!dir.exists(docs_folder)) {
  
  dir.create(docs_folder)
  
}

md_title <- args[4]


# Local debug - do not uncomment
# 
# moba_scores_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/pc_scores.sscore"
# hgdp_scores_file <- "/mnt/archive/snpQc/pc_loadings/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz"
# md_file <- "/mnt/work/marc/github/mobaGenetics-qc/qc-pipeline/docs/snp014/pca_hgdp_moba.md"
# md_title <- "Principal Component Analysys (PCA) in snp014 vs. HDGP"
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
hgdp_scores <- read.table(
  file = hgdp_scores_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  rename(
    iid = s
  ) %>% 
  mutate(
    cohort = "HGDP"
  )

moba_scores <- read.table(
  file = moba_scores_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  mutate(
    cohort = "MoBa"
  )

for (pc_i in 1:20) {
  
  col_name <- glue("pc{pc_i}_sum")
  moba_scores[[col_name]] <- moba_scores[[col_name]] / sqrt(moba_scores$allele_ct)
  
}

plot_data <- rbind(
  hgdp_scores,
  moba_scores %>% 
    select(
      iid,
      pc1 = pc1_sum,
      pc2 = pc2_sum,
      pc3 = pc3_sum,
      pc4 = pc4_sum,
      pc5 = pc5_sum,
      pc6 = pc6_sum,
      pc7 = pc7_sum,
      pc8 = pc8_sum,
      pc9 = pc9_sum,
      pc10 = pc10_sum,
      pc11 = pc11_sum,
      pc12 = pc12_sum,
      pc13 = pc13_sum,
      pc14 = pc14_sum,
      pc15 = pc15_sum,
      pc16 = pc16_sum,
      pc17 = pc17_sum,
      pc18 = pc18_sum,
      pc19 = pc19_sum,
      pc20 = pc20_sum,
      cohort
    )
)

# Merge

cohort_order <- c("HGDP", "MoBa")

plot_data <- plot_data %>% 
  mutate(
    cohort_factor = factor(cohort, levels = cohort_order)
  ) %>% 
  arrange(
    desc(cohort_factor)
  )


# Write docs

write(
  x = paste0("# ", md_title),
  file = md_file,
  append = F
)

write(
  x = paste0("Principal component analysis of the MoBa samples projected over PCs from the Human Genome Diversity Project (HGDP). Based on the analysis protocol of the _COVID-19 Host Genetics Initiative_, see [flagship paper](doi.org/10.1101/2021.03.10.21252820), with adaptations."),
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

n_samples <- as.data.frame(table(plot_data$cohort_factor))

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
  n = length(cohort_order) - 1,
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
      cohort == "MoBa"
    )
  
  kg_data <- plot_data %>% 
    filter(
      cohort != "MoBa"
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
        col = cohort_factor
      )
    ) +
    geom_xsidedensity(
      data = plot_data,
      mapping = aes(
        x = x,
        y = after_stat(density),
        fill = cohort_factor
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = plot_data,
      mapping = aes(
        x = after_stat(density),
        y = y,
        fill = cohort_factor
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
  
  file_name <- paste0(pc_name_x, "_", pc_name_y, "_hgdp.png")
  
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
