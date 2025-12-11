
library(janitor)
library(tidyr)
library(dplyr)
library(ggplot2)
debug <- T
if (debug){
args <- c("/mnt/archive2/moba_genotypes_resources/HRC/",
        "/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod6-imputation/snp007/",
        "/mnt/work/oystein/tmp/imputation_report/imputation_report.md",
        "snp007"
        )
} else {
    args <- commandArgs(TRUE)
}

ref_folder <- args[1]
imputation_folder <- args[2]
md <- args[3]
md_folder <- dirname(md)
if(!dir.exists(md_folder)){
    dir.create(md_folder, recursive = T)
}

write(
  x = paste0("#Imputation report, batch ", args[4]),
  file = md,
  append = F
)


for (chr in c(1:22, "X")) { 
    write(
    x = paste0("## Chromosome ", chr),
    file = md,
    append = T
    )
    ref_file <- paste0(ref_folder, "HRC_af_chr", chr)
    imp_file <- paste0(imputation_folder, "mod6_impute.chr", chr, ".imputed.vcf.gz.info")
    if(!file.exists(ref_file) | !file.exists(imp_file)){
        next
    }
    ref <- read.table(
      file = ref_file,
      header = F,
      col.names = c("chr", "pos", "id", "af"),
      stringsAsFactors = F
    )


    imp <- read.table(
      file = imp_file,
      header = F,
        col.names = c("chr", "pos", "id", "imp", "ref", "alt", "dr2", "af"),
      stringsAsFactors = F
    ) 

    #batr plot of dr2 values
    ggplot(imp, aes(x=af, y=dr2)) + 
        geom_point(alpha=0.3) +
        xlab("Allele frequency") +
        ylab("Dosage R^2") +
        ggtitle(paste0("Chromosome ", chr))
    ggsave(filename = paste0(md_folder, "/imputation_dr2_chr", chr, ".png"), dpi=200, width = 5, height = 5)
    write(
    x = paste0("![](imputation_dr2_chr", chr, ".png)"),
    file = md,
    append = T
    )
    af <- ref %>% select(chr, pos, af_ref = af) %>% inner_join(imp %>% select(chr, pos, af_imp = af), by = c("chr", "pos"))


    ggplot(af, aes(x=af_ref, y=af_imp)) + 
        geom_point(alpha=0.3) +
        geom_abline(slope=1, intercept=0, color="red") +
        xlab("Allele frequency in reference panel") +
        ylab("Allele frequency in imputed data") +
        ggtitle(paste0("Chromosome ", chr))
    ggsave(filename = paste0(md_folder, "/imputation_af_chr", chr, ".png"), dpi=300, width = 5, height = 5)
    write(
    x = paste0("![](imputation_af_chr", chr, ".png)"),
    file = md,
    append = T
    )
}