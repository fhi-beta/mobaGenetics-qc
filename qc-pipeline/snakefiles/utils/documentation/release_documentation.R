
##
#
# This script writes a redme for the release.
#
##

set.seed(11111)


# Libraries

library(glue)


# Command line arguments

args <- commandArgs(TRUE)

release_number <- args[1]

docs_folder <- args[2]

if (!dir.exists(docs_folder)) {
  
  stop(glue("Docs folder {docs_folder} not found."))
  
}

md_name <- args[3]


# Write readme

md_file <- file.path(docs_folder, md_name)

write(
  x = "## Moba genotypes QC",
  file = md_file,
  append = F
)

write(
  x = glue("Documentation of the quality control of the genotypes from the [Norwegian Mother, Father and Child Cohort Study (MoBa)](https://www.fhi.no/en/studies/moba/), release {release_number}."),
  file = md_file,
  append = T
)

write(
  x = "### Post-imputation",
  file = md_file,
  append = T
)

write(
  x = "Quality control conducted post-imputation on the entire release",
  file = md_file,
  append = T
)

write(
  x = "- [Relatedness and sex assignment](mod8_psam_reconstruction.md): Summary statistics on inferred familial relationship and sex assignment. This was conducted using a set of approximately 500,000 high quality genotyped or imputed autosome markers pruned for LD.",
  file = md_file,
  append = T
)
write(
  x = "- [Admixture](pca_1kg_moba.md): Summary statistics on population clustering after running a principal component analysis on the genotypes merged with the 1,000 genomes. This was conducted using a set of approximately 500,000 high quality genotyped or imputedence autosome markers pruned for LD.",
  file = md_file,
  append = T
)

write(
  x = "### Pre-imputation",
  file = md_file,
  append = T
)

write(
  x = "Quality control conducted pre-imputation on individual genotyping batches",
  file = md_file,
  append = T
)

for (batch in list.files(docs_folder)) {
  
  if (startsWith(batch, "snp")) {
    
    batch_folder = file.path(docs_folder, batch)
    
    if (dir.exists(batch_folder)) {
      
      write(
        x = glue("#### {batch}"),
        file = md_file,
        append = T
      )
      
      write(
        x = glue("- [Relatedness and sex assignment]({batch}/fam_reconstruction.md): Summary statistics on inferred familial relationship and sex assignment. This was conducted using a set of high quality autosome markers pruned for LD."),
        file = md_file,
        append = T
      )
      write(
        x = "- [Admixture]({batch}/pca_1kg_moba.md): Summary statistics on population clustering after running a principal component analysis on the genotypes merged with the 1,000 genomes. This was conducted using a set of approximately 500,000 high quality genotyped or imputedence autosome markers pruned for LD.",
        file = md_file,
        append = T
      )
      
    }
    
  }
  
}






