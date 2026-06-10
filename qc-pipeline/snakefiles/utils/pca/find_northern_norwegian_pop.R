library(dplyr)
library(janitor)
pcs_file <- "/mnt/archive3/snpQc/pipeOut_dev/2025.09.25/mod8-release_annotation/mod8_pca_both.pcs"

pcs <- read.table(
  file = pcs_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
) %>% 
  clean_names()


nn <- pcs %>% arrange(pc3) %>% slice_head(n=500) %>% select(iid)


write.table(
  x = nn,
  file = "/mnt/archive2/moba_genotypes_resources/phenotypes/northern_norwegians",
  quote = F,
  row.names = F,
  col.names = F
)