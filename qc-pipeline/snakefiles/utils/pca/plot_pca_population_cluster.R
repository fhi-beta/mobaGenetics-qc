
##
#
# This script plots the results of a PCA and infer the cluster inference of MoBa samples
#
##

# Command line arguments

args <- commandArgs(TRUE)

if (length(args) != 5) {
  
  stop(paste0("Four arguments expected: pcs file, thousand genomes population file, docs folder, md title, output file. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

pcs_file <- args[1]

if (!file.exists(pcs_file)) {
  
  stop("PCs file not found")
  
}

thousand_genomes_populations_file <- args[2]

if (!file.exists(thousand_genomes_populations_file)) {
  
  stop("Thousand genomes population file not found")
  
}

docs_folder <- args[3]

if (!dir.exists(docs_folder)) {
  
  dir.create(docs_folder)
  
}

md_title <- args[4]
export_file <- args[5]


# Libraries
library(janitor)
library(tidyr)
library(dplyr)
library(ggplot2)
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

md_file <- file.path(docs_folder, "pca.md")

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

n_samples <- as.data.frame(table(plot_data$pop))

for (sample_i in 1:nrow(n_samples)) {
  
  write(
    x = paste0("| ", n_samples[sample_i, 1], " | ", n_samples[sample_i, 2], " |"),
    file = md_file,
    append = T
  )
  
}


# Plot the PCs against each other

for (pc_i in 1:9) {
  
  pc_name_x <- paste0("pc", pc_i)
  pc_name_y <- paste0("pc", pc_i + 1)
  
  plot_data$x <- plot_data[[pc_name_x]]
  plot_data$y <- plot_data[[pc_name_y]]
  
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
      data = plot_data,
      mapping = aes(
        x = x,
        y = y,
        col = pop
      ),
      alpha = 0.5
    ) +
    scale_x_continuous(
      name = pc_name_x
    ) +
    scale_y_continuous(
      name = pc_name_y
    )
  
  png(
    filename = paste0("qc-pipeline/docs/snp012/pca/", pc_name_x, "_", pc_name_y, ".png"),
    width = 800,
    height = 600
  )
  grid.draw(plot)
  device <- dev.off()
  
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
