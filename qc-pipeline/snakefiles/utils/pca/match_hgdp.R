
##
#
# This script matches the variants in a batch to the loadings of the HGDP and 1k genomes project.
#
##

# Command line arguments

args <- commandArgs(TRUE)

if (length(args) != 4) {
  
  stop(paste0("Four arguments expected: pcs file, thousand genomes population file, md file, md title, output file. ", length(args), " found: ", paste(args, collapse = ", ")))
  
}

variant_file <- args[1]

if (!file.exists(variant_file)) {
  
  stop("Variant file not found")
  
}

loadings_file <- args[2]

if (!file.exists(loadings_file)) {
  
  stop("PC loadings not found")
  
}

proxy_cache_folder <- args[3]

export_file <- args[4]


# DEBUG
# variant_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/mod3_convert_plink2.pvar"
# loadings_file <- "/mnt/archive/snpQc/pc_loadings/hgdp_tgp_pca_covid19hgi_snps_loadings.GRCh37.plink.tsv"
# proxy_cache_folder <- "/mnt/archive/snpQc/pipeOut/tmp/ldlink_cache"
# export_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/loadings_hdpg_1kg"


# Libraries

library(janitor)
library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(LDlinkR)


# Load the data

variant_table <- read.table(
  file = variant_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ""
) %>% 
  clean_names() %>% 
  rename(
    chr = x_chrom
  )

loadings_table <- read.table(
  file = loadings_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

original_names <- names(loadings_table)

loadings_table <- loadings_table %>% 
  clean_names()


# Match the variants

matched_loadings <- list()

for (variant_i in 1:nrow(variant_table)) {
  
  variant_id <- variant_table$id[variant_i]
  variant_chr <- variant_table$chr[variant_i]
  variant_pos <- variant_table$pos[variant_i]
  variant_ref <- variant_table$ref[variant_i]
  variant_alt <- variant_table$alt[variant_i]
  
  print(glue("Processing {variant_id} ({variant_i} of {nrow(variant_table)})"))
  
  temp_id <- paste(variant_chr, variant_pos, variant_ref, variant_alt, sep = ":")
  
  loading_i <- which(loadings_table$id == temp_id)
  
  if (length(loading_i) > 0) {
    
    if (length(loading_i) > 1) {
      
      stop("Duplicate snp")
      
    }
    
    temp <- loadings_table[loading_i, ]
    temp$id <- variant_id
    
    matched_loadings[[length(matched_loadings) + 1]] <- temp
    temp$alt <- variant_alt
    
  } else {
    
    temp_id <- paste(variant_chr, variant_pos, variant_alt, variant_ref, sep = ":")
    
    loading_i <- which(loadings_table$id == temp_id)
    
    if (length(loading_i) > 0) {
      
      if (length(loading_i) > 1) {
        
        stop("Duplicate snp")
        
      }
      
      variant_id <- variant_table$id[variant_i]
      temp <- loadings_table[loading_i, ]
      temp$id <- variant_id
      temp$alt <- variant_ref
      
      matched_loadings[[length(matched_loadings) + 1]] <- temp
      
    } else {
      
      id_is_rsid <- startsWith(variant_id, "rs") && is.finite(as.numeric(substring(variant_id, 3)))
      
      if (id_is_rsid) {
        
        temp_id <- variant_id
        cache_file_name <- paste0(variant_id, ".gz")
        
      } else {
        
        temp_id <- paste0("chr", variant_chr, ":", variant_pos)
        cache_file_name <- paste0("chr", variant_chr, "_", variant_pos, ".gz")
        
      }
      
      cache_file <- file.path(proxy_cache_folder, cache_file_name)
      
      if (!file.exists(cache_file)) {
        
        proxy_table <- LDproxy(
          snp = temp_id, 
          pop = "ALL", 
          token = "972f33fe5966"
        ) %>% 
          clean_names()  %>% 
          filter(
            r2 >= 0.2
          ) %>% 
          select(
            rs_number,
            coord,
            alleles,
            distance,
            dprime,
            r2,
            correlated_alleles
          )
        
        write.table(
          x = proxy_table,
          file = gzfile(cache_file),
          sep = "\t",
          col.names = T,
          row.names = F
        )
        
      } else {
        
        proxy_table <- read.table(
          file = cache_file,
          header = T,
          sep = "\t",
          stringsAsFactors = F
        )
        
      }
      
      if (nrow(proxy_table) > 0) {
        
        proxy_table <- proxy_table %>% 
          arrange(desc(r2), desc(abs(distance)))
        
        for (proxy_i in 1:nrow(proxy_table)) {
          
          correlated_alleles <- proxy_table$correlated_alleles[proxy_i]
          
          temp_alleles <- strsplit(correlated_alleles, ",")[[1]]
          temp_alleles_1 <- strsplit(temp_alleles[1], "=")[[1]]
          temp_alleles_2 <- strsplit(temp_alleles[2], "=")[[1]]
          temp_alleles_ref <- temp_alleles_1[1]
          proxy_ref <- temp_alleles_1[2]
          temp_alleles_alt <- temp_alleles_2[1]
          proxy_alt <- temp_alleles_2[2]
          
          if (temp_alleles_ref == variant_ref && temp_alleles_alt == variant_alt ||
              temp_alleles_ref == variant_alt && temp_alleles_alt == variant_ref) {
            
            aligned_alleles_1 <- paste0(variant_ref, "=", proxy_ref, ",", variant_alt, "=", proxy_alt)
            aligned_alleles_2 <- paste0(variant_alt, "=", proxy_alt, ",", variant_ref, "=", proxy_ref)
            
            swapped_alleles_1 <- paste0(variant_ref, "=", proxy_alt, ",", variant_alt, "=", proxy_ref)
            swapped_alleles_2 <- paste0(variant_alt, "=", proxy_ref, ",", variant_ref, "=", proxy_alt)
            
            if (swapped_alleles_1 == correlated_alleles || swapped_alleles_2 == correlated_alleles) {
              
              proxy_ref <- temp_alleles_2[2]
              proxy_alt <- temp_alleles_1[2]
              
            } else if (aligned_alleles_1 != correlated_alleles && aligned_alleles_2 != correlated_alleles) {
              
              stop("Allele mismatch")
              
            }
            
            proxy_id <- str_replace(proxy_table$coord[proxy_i], "chr", "")
            temp_id <- paste(proxy_id, proxy_ref, proxy_alt, sep = ":")
            
            loading_i <- which(loadings_table$id == temp_id)
            
            if (length(loading_i) > 0) {
              
              if (length(loading_i) > 1) {
                
                stop("Duplicate snp")
                
              }
              
              temp <- loadings_table[loading_i, ]
              temp$id <- variant_id
              
              matched_loadings[[length(matched_loadings) + 1]] <- temp
              temp$alt <- variant_alt
              
            } else {
              
              temp_id <- paste(proxy_id, proxy_alt, proxy_ref, sep = ":")
              
              loading_i <- which(loadings_table$id == temp_id)
              
              if (length(loading_i) > 0) {
                
                if (length(loading_i) > 1) {
                  
                  stop("Duplicate snp")
                  
                }
                
                temp <- loadings_table[loading_i, ]
                temp$id <- variant_id
                
                matched_loadings[[length(matched_loadings) + 1]] <- temp
                temp$alt <- variant_ref
                
              }
            }
          }
        }
      }
    }
  }
}

matched_loadings <- do.call("rbind", matched_loadings)

print(glue("{nrow(matched_loadings)} variants matched to loading ({nrow(variant_table)} variants in MoBa, {nrow(loadings_table)} variants in loadings)"))


# Write results

names(matched_loadings) <- original_names

write.table(
  x = matched_loadings,
  file = export_file,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

  