

debug <- T
if (debug) {
  args <- c("/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_pedigree_ibd_estimate.kin0",
  "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_batch_table_batch")
} else {
  args <- commandArgs(TRUE)
}


library(igraph)
library(dplyr)
library(ggplot2)
library(grid)

genomic_relatedness_table <- read.table(
  file = args[1],
  header = T,
  stringsAsFactors = F
)

batch_table <- read.table(
  file = args[2],
  header = T,
  stringsAsFactors = F
)

sample_ids <- batch_table$iid

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
  family = paste0("kinship_cluster_", nuclear_families$membership),
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
id_to_family <- id_to_family %>% select(fid = family, iid = id)
families <- unique(id_to_family$fid)
for (f in families) {
  family_members <- subset(id_to_family, fid == f) %>% select(iid)
  family_file <- paste0("/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/families/", f, ".psam")
  write.table(
    family_members,
    file = family_file,
    quote = F,
    row.names = F,
    col.names = F
  )
}