##
#
# This script formats the results from regenie for further analysis using LocusZoom and COJO.
# 
##


print(paste0(Sys.time(), " - Formatting regenie results for further analysis"))


# Command line input

args <- commandArgs(TRUE)

if (length(args) != 7) {
  
  stop(paste0("Seven parameters expected. ", length(args), " found."))
  
}

regenie_file <- args[1]

if (!file.exists(regenie_file)) {
  
  stop(paste0("Regenie file ", regenie_file, " not found."))
  
}

regenieFileName <- tools::file_path_sans_ext(basename(regenie_file))

output_folder <- args[2]

if (!dir.exists(output_folder)) {
  
  stop(paste0("Output folder ", output_folder, " not found."))
  
}

done_file <- args[3]

pheno_file <- args[4]

if (!file.exists(pheno_file)) {
  
  stop(paste0("Pheno file ", regenie_file, " not found."))
  
}

pheno_name <- args[5]

lz_folder <- args[6]

if (!dir.exists(lz_folder)) {
  
  stop(paste0("Locus zoom folder ", lz_folder, " not found."))
  
}

threshold <- args[7]

if (is.na(as.numeric(threshold)))  {

  stop(paste0("Threshold ", threshold, " could not be parsed as a number."))

}

threshold <- as.numeric(threshold)


# Read phenotypes

print(paste0(Sys.time(), " - Formatting ", regenieFileName, ": reading phenotypes"))

phenotypes <- read.table(
  file = pheno_file,
  header = T,
  sep = " ",
  stringsAsFactors = F
)

if (!pheno_name %in% names(phenotypes)) {
  
  stop(paste0(pheno_name, " not found in ", pheno_file))
  
}

n_samples <- sum(!is.na(phenotypes[[pheno_name]]))


# Read results from regenie

print(paste0(Sys.time(), " - Formatting ", regenieFileName, ": reading association results"))

regenie_results <- read.table(
  file = regenie_file,
  header = T,
  sep = " ",
  stringsAsFactors = F
)

regenie_results <- regenie_results[is.finite(regenie_results$LOG10P), ]


# Export in a LocusZoom-friendly format

print(paste0(Sys.time(), " - Formatting ", regenieFileName, ": exporting for LocusZoom"))

lz_export <- data.frame(
  id = regenie_results$ID,
  chromosome = regenie_results$CHROM,
  position = regenie_results$GENPOS,
  tested_allele = regenie_results$ALLELE1,
  other_allele = regenie_results$ALLELE0,
  tested_allele_frequency = regenie_results$A1FREQ,
  beta = regenie_results$BETA,
  se = regenie_results$SE,
  p = 10^-regenie_results$LOG10P,
  stringsAsFactors = F
)

output_file <- file.path(lz_folder, paste0(regenieFileName, ".lz.gz"))

write.table(
  x = lz_export,
  file = gzfile(output_file),
  col.names = T,
  row.names = F,
  quote = F,
  sep = "\t"
)



# Write one cojo input file per chromosome

chr_22 <- F

for (chromosome in unique(regenie_results$CHROM)) {
  
  if (chromosome == 22) {
    
    chr_22 <- T
    
  }
  
  print(paste0(Sys.time(), " - Formatting ", regenieFileName, ": exporting results for COJO chromosome ", chromosome))
  
  output_file <- file.path(output_folder, paste0(regenieFileName, ".chr", chromosome, ".ma"))
  
  chromosome_results <- regenie_results[regenie_results$CHROM == chromosome, ]
  
  chr_export <- data.frame(
    SNP = chromosome_results$ID,
    A1 = chromosome_results$ALLELE1,
    A2 = chromosome_results$ALLELE0,
    freq = chromosome_results$A1FREQ,
    b = chromosome_results$BETA,
    se = chromosome_results$SE,
    p = 10^-chromosome_results$LOG10P,
    N =  n_samples,
    stringsAsFactors = F
  )
  
  if (min(chr_export$p) < 0 | max(chr_export$p) > 1) {
    
    stop(paste0("p-value out of range: min ", min(chr_export$p), ", max ", max(chr_export$p)))
    
  }
  
  write.table(
    x = chr_export,
    file = output_file,
    col.names = T,
    row.names = F,
    quote = F,
    sep = " "
  )
  
  
  output_file <- file.path(output_folder, paste0(regenieFileName, ".chr", chromosome, ".snpList"))
  
  snpList <- chr_export[, c("SNP")]
  
  write.table(
    x = snpList,
    file = output_file,
    col.names = F,
    row.names = F,
    quote = F,
    sep = " "
  )
  
  output_file <- file.path(output_folder, paste0(regenieFileName, ".chr", chromosome, ".no_snp"))
  
  if (nrow(chr_export) == 0 || min(chr_export$p) > threshold) {
    
    writeLines(
      text = "No GWAS-significant SNP", 
      output_file
    )
    
  } else if (file.exists(output_file)) {
    
    file.remove(output_file)
    
  }
}

if (!chr_22) {
  
  stop(paste0("Chromosome 22 not found in ", regenieFileName, " please check that the regenie output file is not trucated."))
  
}


# Tell the pipeline we are done

writeLines(
  text = paste0(Sys.time(), ": Done"), 
  done_file
)



