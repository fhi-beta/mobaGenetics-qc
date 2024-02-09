
##
#
# This script builds sets of unrelated participants.
# It is based on the `accumPIHAT` script written by Jonas on 2017.01.15.
#
##

set.seed(11111)


# Command line arguments
args <- commandArgs(TRUE)

if (length(args) != 8) {
  
  stop(paste0("Eight arguments expected: genome file, sex check file, registry file, current fam file, destination file, exclusion file, md file, title. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

king_file <- args[1]

if (!file.exists(king_file)) {
  
  stop("King file not found")
  
}

kinship_threshold <- args[2]

if (is.na(as.numeric(kinship_threshold))) {
  
  stop(psate0("Input for `kinship_threshold` (", kinship_threshold, ") could not be parsed as a number."))
  
}

accumulated_kinship_threshold <- args[3]

if (is.na(as.numeric(accumulated_kinship_threshold))) {
  
  stop(psate0("Input for `accumulated_kinship_threshold` (", accumulated_kinship_threshold, ") could not be parsed as a number."))
  
}

exclusion_file <- args[4]

md_file <- args[5]

title <- args[6]


# Libraries

library(igraph)
library(dplyr)
library(ggplot2)
library(grid)

theme_set(theme_bw(base_size = 24))


# Load data

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
  x = "## Relatedness",
  file = md_file,
  append = T
)

ibd_plot <- ggplot() +
  geom_point(
    data = genomic_relatedness_table,
    mapping = aes(
      x = IBS0,
      y = IBD2Seg,
      col = relationship
    ),
    alpha = 0.8
  ) +
  scale_x_continuous(
    name = "IBS0 [Porportion of SNPs with zero IBS]"
  ) +
  scale_y_continuous(
    name = "Kinship"
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


# Export ids of samples to exclude

to_remove_fam <- fam_data

write.table(
  x = to_remove_fam,
  file = exclusion_file,
  append = F,
  col.names = F,
  row.names = F,
  sep = " ",
  quote = F
)

