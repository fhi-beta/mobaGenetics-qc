
##
#
# This script matches the variants in a batch to the loadings of the HGDP and 1k genomes project.
#
##

# Command line arguments

args <- commandArgs(TRUE)


if (length(args) != 6) {

  stop(paste0("Four arguments expected: (1) variant file, (2) PC loadings file, (3) proxies cache file, (4) proxies database, (5) file where to export the matched loadings, (6) file where to export the ids of the matched variants. ", length(args), " found: ", paste(args, collapse = ", ")))

}

variant_file <- args[1]

if (!file.exists(variant_file)) {

  stop("Variant file not found")

}

loadings_file <- args[2]

if (!file.exists(loadings_file)) {

  stop("PC loadings not found")

}

frequency_file <- args[3]

if (!file.exists(frequency_file)) {
  
  stop("PC loadings not found")
  
}

proxies_cache_stem <- args[4]

proxy_db <- args[5]

loading_export_file <- args[6]

frequencies_export_file <- args[7]

variants_export_file <- args[8]


# DEBUG
# variant_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/mod3_convert_plink2.pvar"
# loadings_file <- "/mnt/archive/snpQc/pc_loadings/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.tsv"
# frequency_file <- "/mnt/archive/snpQc/pc_loadings/hgdp_tgp_pca_covid19hgi_snps_loadings.rsid.plink.afreq"
# proxies_cache_stem <- "/mnt/archive/snpQc/pc_loadings/proxies_cache"
# proxy_db <- "/mnt/archive/topld/db/ld_db"
# loading_export_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/loadings_hdpg_1kg_proxies"
# frequencies_export_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/frequencies_hdpg_1kg_proxies"
# variants_export_file <- "/mnt/archive/snpQc/pipeOut_dev/mod3-good-markers/snp014/variants_hdpg_1kg"


# Libraries

library(janitor)
library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(DBI)
library(dbplyr)


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

loadings_header <- names(loadings_table)

loadings_table <- loadings_table %>% 
  clean_names()

frequency_table <- read.table(
  file = frequency_file,
  header = T,
  sep = "\t",
  stringsAsFactors = F,
  comment.char = ""
)
names(frequency_table) <- c("id", "ref", "alt", "freq")


# Set up database and caches for LD

db_connection <- dbConnect(RSQLite::SQLite(), proxy_db)

annotation_table <- NULL
proxy_tables <- NULL

proxies_cache <- data.frame(
  moba_snp = character(0),
  moba_ref = character(0),
  moba_alt = character(0),
  loading_proxy = character(0),
  allele_swap = logical(0),
  stringsAsFactors = F
)

no_proxies_cache <- NULL

superpopulations <- c("AFR", "EAS", "EUR", "SAS")


# Match the variants
# Note: Currently MoBa is GRCh37, TopLD is GRCh38, to avoid lifts we map by rsid only. Once MoBa is GRCh38 we can map by variant coordinates and alleles.

matched_loadings <- list()
matched_frequencies <- list()
new_proxies <- list()

current_chr_proxies <- -1
current_chr_cache <- -1

for (variant_i in 1:nrow(variant_table)) {
  
  variant_id <- variant_table$id[variant_i]
  variant_chr <- variant_table$chr[variant_i]
  variant_pos <- variant_table$pos[variant_i]
  variant_ref <- variant_table$ref[variant_i]
  variant_alt <- variant_table$alt[variant_i]
  
  if (variant_chr != current_chr_cache) {
    
    if (nrow(proxies_cache) > 0) {
      
      proxies_cache <- proxies_cache %>% distinct()
      
      proxies_cache_file <- paste0(proxies_cache_stem, "_", current_chr_cache, ".gz")
      
      write.table(
        x = proxies_cache,
        file = gzfile(proxies_cache_file),
        sep = "\t",
        col.names = T,
        row.names = F,
        quote = F
      )
      
    }
    
    if (length(no_proxies_cache) > 0) {
      
      no_proxies_cache <- unique(no_proxies_cache)
      
      no_proxies_cache_file <- paste0(proxies_cache_stem, "_no_proxies_", current_chr_cache, ".gz")
      
      writeLines(no_proxies_cache, con = gzfile(no_proxies_cache_file))
      
    }
    
    proxies_cache_file <- paste0(proxies_cache_stem, "_", variant_chr, ".gz")
    
    if (file.exists(proxies_cache_file)) {
      
      proxies_cache <- read.table(
        file = proxies_cache_file,
        header = T,
        sep = "\t",
        stringsAsFactors = F
      )
      
    } else {
      
      proxies_cache <- data.frame(
        moba_snp = character(0),
        moba_ref = character(0),
        moba_alt = character(0),
        loading_proxy = character(0),
        allele_swap = logical(0),
        stringsAsFactors = F
      )
      
    }
    
    no_proxies_cache_file <- paste0(proxies_cache_stem, "_no_proxies_", variant_chr, ".gz")
    
    if (file.exists(no_proxies_cache_file)) {
      
      no_proxies_cache <- readLines(no_proxies_cache_file)
      
    } else {
      
      no_proxies_cache <- c()
      
    }
    
    current_chr_cache <- variant_chr
    
  }
  
  # print(glue("{Sys.time()}    Processing {variant_id} ({variant_i} of {nrow(variant_table)})"))
  
  rs_id <- variant_id
  
  if (startsWith(rs_id, "GSA-")) {
    
    rs_id <- substring(rs_id, 5)
    
  }
  
  if (startsWith(rs_id, "BOT2-")) {
    
    rs_id <- substring(rs_id, 6)
    
  }
  
  if (startsWith(rs_id, "seq-")) {
    
    rs_id <- substring(rs_id, 5)
    
  }
  
  if (startsWith(rs_id, "DUP-")) {
    
    rs_id <- substring(rs_id, 5)
    
  }
  
  if (startsWith(rs_id, "TOP-")) {
    
    rs_id <- substring(rs_id, 5)
    
  }
  
  if (startsWith(rs_id, "seq_")) {
    
    rs_id <- substring(rs_id, 5)
    
  }
  
  if (endsWith(rs_id, "_dup")) {
    
    rs_id <- substring(rs_id, 1, nchar(rs_id) - 4)
    
  }
  
  if (!startsWith(x = rs_id, prefix = "rs") & grepl(x = rs_id, pattern = "rs", fixed = TRUE)) {
    
    stop(glue("Variant '{rs_id}' not supported."))
    
  }
  
  if (startsWith(x = rs_id, prefix = "rs") && !rs_id %in% no_proxies_cache) {
    
    loading_i <- which(loadings_table$id == rs_id)
    
    if (length(loading_i) > 0) {
      
      if (length(loading_i) > 1) {
        
        stop("Duplicate snp")
        
      }
      
      temp <- loadings_table[loading_i, ]
      temp$id[1] <- variant_id
      temp$alt[1] <- variant_alt
      
      matched_loadings[[length(matched_loadings) + 1]] <- temp
      
      frequency_i <- which(frequency_table$id == rs_id)
      
      if (length(loading_i) != 1) {
        
        stop("Frequency mismatch")
        
      }
      
      temp <- frequency_table[loading_i, ]
      temp$id[1] <- variant_id
      temp$ref[1] <- variant_ref
      temp$alt[1] <- variant_alt
      
      matched_frequencies[[length(matched_frequencies) + 1]] <- temp
      
    } else if (rs_id %in% proxies_cache$moba_snp) {
      
      proxy_i <- which(proxies_cache$moba_snp == rs_id)[1]
      
      allele_swap <- proxies_cache$allele_swap[proxy_i]
      
      if (variant_alt == proxies_cache$moba_ref[proxy_i] && variant_ref == proxies_cache$moba_alt[proxy_i]) {
        
        allele_swap <- !allele_swap
        
      } else if (variant_alt != proxies_cache$moba_alt[proxy_i] || variant_alt == proxies_cache$moba_ref[proxy_i]) {
        
        stop("Allele mismatch between MoBa versions")
        
      }
      
      proxy_rs_id <- proxies_cache$loading_proxy[proxy_i]
      
      loading_i <- which(loadings_table$id == proxy_rs_id)
      
      if (length(loading_i) != 1) {
        
        stop("Incorrect number of loadings")
        
      }
      
      temp <- loadings_table[loading_i, ]
      temp$id[1] <- variant_id
      
      if (!allele_swap) {
        
        temp$alt[1] <- variant_alt
        
      } else {
        
        temp$alt[1] <- variant_ref
        
      }
      
      matched_loadings[[length(matched_loadings) + 1]] <- temp
      
      frequency_i <- which(frequency_table$id == proxy_rs_id)
      
      if (length(frequency_i) != 1) {
        
        stop("Incorrect number of frequencies")
        
      }
      
      temp <- frequency_table[frequency_i, ]
      temp$id[1] <- variant_id
        temp$ref[1] <- variant_ref
        temp$alt[1] <- variant_alt
      
      if (allele_swap) {
        
        temp$freq[1] <- 1.0 - temp$freq[1]
        
      }
      
      matched_frequencies[[length(matched_frequencies) + 1]] <- temp
      
    } else {
      
      if (variant_chr != current_chr_proxies) {
        
        print(glue("{Sys.time()}    Loading LD annotation for chromosome {variant_chr}"))
        
        annotation_table <- dbReadTable(db_connection, glue("annotation_{variant_chr}"))
        
        proxy_tables <- list()
        
        for (superpopulation in superpopulations) {
          
          print(glue("{Sys.time()}    Loading LD values for chromosome {variant_chr} superpopulation {superpopulation}"))
          
          proxy_tables[[superpopulation]] <- dbReadTable(db_connection, glue("ld_{superpopulation}_{variant_chr}"))
          
        }
        
        current_chr_proxies <- variant_chr
        
      }
      
      proxy_found <- F
      
      annotation_is <- which(annotation_table$rs_id == rs_id)
      
      if (length(annotation_is) > 0) {
        
        for(annotation_i in annotation_is) {
          
          topld_id <- annotation_table$uniq_id[annotation_i]
          topld_ref <- annotation_table$ref[annotation_i]
          topld_alt <- annotation_table$alt[annotation_i]
          
          proxies <- list()
          
          for (superpopulation in superpopulations) {
          
          population_proxies <- proxy_tables[[superpopulation]][proxy_tables[[superpopulation]]$uniq_id_1 == topld_id | proxy_tables[[superpopulation]]$uniq_id_2 == topld_id, ]
          
          if (nrow(population_proxies) > 0) {
            
            population_proxies$superpopulation <- superpopulation
            
          } else {
            
            population_proxies <- data.frame(
              snp1 = c(),
              snp2 = c(),
              uniq_id_1 = c(),
              uniq_id_2 = c(),
              r2 = c(),
              dprime = c(),
              x_corr = c(),
              superpopulation = c()
            )
          }
          
          proxies[[superpopulation]] <- population_proxies
          
          }
          
          proxies <- do.call(rbind, proxies)
          
          if (nrow(proxies) > 0) {
            
            proxies <- proxies %>% 
              arrange(
                desc(r2)
              )
            
            for (proxy_i in 1:nrow(proxies)) {
              
              if (proxies$uniq_id_1[proxy_i] == topld_id) {
                
                proxy <- proxies$uniq_id_2[proxy_i]
                
              } else {
                
                proxy <- proxies$uniq_id_1[proxy_i]
                
              }
              
              proxy_annotation_is <- which(annotation_table$uniq_id == proxy)
              
              for (proxy_annotation_i in proxy_annotation_is) {
                
                proxy_rs_id <- annotation_table$rs_id[proxy_annotation_i]
                
                if (startsWith(proxy_rs_id, "rs")) {
                  
                  loading_i <- which(loadings_table$id == proxy_rs_id)
                  
                  if (length(loading_i) > 0) {
                    
                    if (length(loading_i) > 1) {
                      
                      stop("Duplicate snp")
                      
                    }
                    
                    allele_swap <- proxies$x_corr[proxy_i] != '+'
                    
                    correct_alleles <- F # Check for multi-allelic variants
                    
                    if (variant_alt == topld_ref && variant_ref == topld_alt) {
                      
                      allele_swap <- !allele_swap
                      correct_alleles <- T
                      
                    } else if (variant_alt == topld_alt && variant_ref == topld_ref) {
                      
                      correct_alleles <- T
                      
                    }
                    
                    if (correct_alleles) {
                      
                      temp <- loadings_table[loading_i, ]
                      temp$id[1] <- variant_id
                      
                      if (!allele_swap) {
                        
                        temp$alt[1] <- variant_alt
                        
                      } else {
                        
                        temp$alt[1] <- variant_ref
                        
                      }
                      
                      matched_loadings[[length(matched_loadings) + 1]] <- temp
                      
                      frequency_i <- which(frequency_table$id == proxy_rs_id)
                      
                      if (length(frequency_i) != 1) {
                        
                        stop("Incorrect number of frequency")
                        
                      }
                      
                      temp <- frequency_table[frequency_i, ]
                      temp$id[1] <- variant_id
                      temp$ref[1] <- variant_ref
                      temp$alt[1] <- variant_alt
                      
                      if (allele_swap) {
                        
                        temp$freq[1] <- 1.0 - temp$freq[1]
                        
                      }
                      
                      matched_frequencies[[length(matched_frequencies) + 1]] <- temp
                      
                      new_proxy <- data.frame(
                        moba_snp = rs_id,
                        moba_ref = variant_ref,
                        moba_alt = variant_alt,
                        loading_proxy = proxy_rs_id,
                        allele_swap = allele_swap,
                        stringsAsFactors = F
                      )
                      
                      new_proxies[[length(new_proxies) + 1]] <- new_proxy
                      
                      proxy_found <- T
                      
                      break
                      
                    }
                  }
                }
              }
            }
          }
        }
      }
      
      if (!proxy_found) {
        
        no_proxies_cache[length(no_proxies_cache) + 1] <- rs_id
        
      } else {
        
        new_proxies <- do.call(rbind, new_proxies)
        proxies_cache <- rbind(proxies_cache, new_proxies)
        
        new_proxies <- list()
        
      }
    }
  }
}

new_proxies <- do.call(rbind, new_proxies)
proxies_cache <- rbind(proxies_cache, new_proxies) %>% distinct()

proxies_cache_file <- paste0(proxies_cache_stem, "_", current_chr_cache, ".gz")

write.table(
  x = proxies_cache,
  file = gzfile(proxies_cache_file),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

no_proxies_cache <- unique(no_proxies_cache)

no_proxies_cache_file <- paste0(proxies_cache_stem, "_no_proxies_", current_chr_cache, ".gz")

writeLines(no_proxies_cache, con = gzfile(no_proxies_cache_file))

print(glue("{Sys.time()}    {nrow(matched_loadings)} variants matched to loading ({nrow(variant_table)} variants in MoBa, {nrow(loadings_table)} variants in loadings)"))


# Write results

matched_loadings_df <- do.call(rbind, matched_loadings)
matched_loadings_df <- matched_loadings_df %>% 
  group_by(
    id, alt
  ) %>% 
  summarize(
    PC1 = mean(pc1),
    pc2 = mean(pc2),
    pc3 = mean(pc3),
    pc4 = mean(pc4),
    pc5 = mean(pc5),
    pc6 = mean(pc6),
    pc7 = mean(pc7),
    pc8 = mean(pc8),
    pc9 = mean(pc9),
    pc10 = mean(pc10),
    pc11 = mean(pc11),
    pc12 = mean(pc12),
    pc13 = mean(pc13),
    pc14 = mean(pc14),
    pc15 = mean(pc15),
    pc16 = mean(pc16),
    pc17 = mean(pc17),
    pc18 = mean(pc18),
    pc19 = mean(pc19),
    pc20 = mean(pc20),
    .groups = "drop"
  )

names(matched_loadings_df) <- loadings_header

write.table(
  x = matched_loadings_df,
  file = loading_export_file,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

matched_frequencies <- do.call(rbind, matched_frequencies)
matched_frequencies <- matched_frequencies %>% 
  group_by(
    id, ref, alt
  ) %>% 
  summarize(
    freq = mean(freq)
  )
names(matched_frequencies) <- c("#ID", "REF", "ALT", "ALT1_FREQ")

write.table(
  x = matched_frequencies,
  file = frequencies_export_file,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

writeLines(
  text = matched_loadings_df$ID,
  con = variants_export_file
)
