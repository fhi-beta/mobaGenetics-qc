
##
#
# This script merges the cojo files from different chromosomes
#
##


print(paste0(Sys.time(), " - Merging COJO results"))


# Command line input

args <- commandArgs(TRUE)

if (length(args) != 2) {
  
   stop(paste0("Two command line argument expected: cojo file path, output. ", length(args), " found."))
  
}

cojo_file_stem <- args[1]
output_file_stem <- args[2]

# Libraries

library(stringr)
library(glue)
library(dplyr)
library(janitor)


# Load data

chr_cma <- list()
chr_jma <- list()

for (chromosome in 1:22) {
  
  cojo_file <- paste0(
    str_replace_all(
    string = cojo_file_stem,
    pattern = "_wildcard_chr",
    replacement = as.character(chromosome)
    ),
    ".cma.cojo"
  )
  
  if (file.exists(cojo_file)) {
    
    chr_data_frame <- read.table(
      file = cojo_file,
      header = T,
      sep = "\t",
      stringsAsFactors = F
    ) %>%
      clean_names() %>% 
      filter(
        p_c < 5e-8
      )
    
    if (nrow(chr_data_frame) > 0) {
      
      chr_cma[[length(chr_cma) + 1]] <- chr_data_frame
      
    }

    file.remove(cojo_file)

  }
  
  cojo_file <- paste0(
    str_replace_all(
      string = cojo_file_stem,
      pattern = "_wildcard_chr",
      replacement = as.character(chromosome)
    ),
    ".jma.cojo"
  )
  
  if (file.exists(cojo_file)) {
    
    chr_data_frame <- read.table(
      file = cojo_file,
      header = T,
      sep = "\t",
      stringsAsFactors = F
    ) %>%
      clean_names()
    
    if (nrow(chr_data_frame) > 0) {
      
      chr_jma[[length(chr_jma) + 1]] <- chr_data_frame
      
    }
  }
}

if (length(chr_cma) > 0) {
  
  cma <- do.call(
    what = "rbind",
    chr_cma
  ) %>% 
    mutate(
      ref_a = ifelse(ref_a == "TRUE", "T", ref_a)
    )
  
} else {
  
  cma <- data.frame(
    chr = character(),
    snp = character(),
    bp = character(),
    ref_a = character(),
    freq = character(),
    b = character(),
    se = character(),
    p = character(),
    n = character(),
    freq_geno = character(),
    b_c = character(),
    b_c_se = character(),
    p_c = character(),
    stringsAsFactors = F
  )
  
}

output_file <- paste0(output_file_stem, ".cma.cojo")

write.table(
  x = cma,
  file = output_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
  )

if (length(chr_jma) > 0) {
  
  jma <- do.call(
    what = "rbind",
    chr_jma
  ) %>% 
    mutate(
      ref_a = ifelse(ref_a == "TRUE", "T", ref_a)
    )
  
} else {
  
  jma <- data.frame(
    chr = character(),
    snp = character(),
    bp = character(),
    ref_a = character(),
    freq = character(),
    b = character(),
    se = character(),
    p = character(),
    n = character(),
    freq_geno = character(),
    b_j = character(),
    b_j_se = character(),
    p_j = character(),
    ld_r = character(),
    stringsAsFactors = F
  )
  
}

output_file <- paste0(output_file_stem, ".jma.cojo")

write.table(
  x = jma,
  file = output_file,
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)


