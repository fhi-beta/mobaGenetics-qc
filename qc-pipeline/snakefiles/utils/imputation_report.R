
library(janitor)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggside)
debug <- F
if (debug){
args <- c("/mnt/archive2/moba_genotypes_resources/HRC",
        "/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod6-imputation/snp008",
        "/mnt/work/oystein/tmp/snp008/imputation_report/imputation_report.md"
        )
} else {
    args <- commandArgs(TRUE)
}

ref_folder <- args[1]
imputation_folder <- args[2]
md <- args[3]
batch <- basename(dirname(md))
md_folder <- dirname(md)
chromosomes <- c(1:22, "X", "PAR1", "PAR2")
if(!dir.exists(md_folder)){
    dir.create(md_folder, recursive = T)
}
plot_folder <- paste0(md_folder, "/imputation_plots/")
if(!dir.exists(plot_folder)){
    dir.create(plot_folder)
}
relative_plot_folder <- "imputation_plots/"

write(
  x = paste0("#Imputation report, batch ", batch),
  file = md,
  append = F
)


for (chr in chromosomes) {
    ref_file <- paste0(ref_folder, "/HRC_af_chr", chr)
    imp_file <- paste0(imputation_folder, "/mod6_impute.chr", chr, ".imputed.vcf.gz.info")
    if(!file.exists(ref_file) | !file.exists(imp_file)){
        if(!file.exists(ref_file)){
            print(paste0("Reference file not found: ", ref_file))
        }
        if(!file.exists(imp_file)){
            print(paste0("Imputation info file not found: ", imp_file))
        }
        next
    }
    write(
    x = paste0("## Chromosome ", chr),
    file = md,
    append = T
    )
    
    ref <- read.table(
      file = ref_file,
      header = F,
      col.names = c("chr", "pos", "id", "af"),
      colClasses = c("character", "numeric", "character", "numeric"),
      stringsAsFactors = F
    )


    imp <- read.table(
      file = imp_file,
      header = F,
      col.names = c("chr", "pos", "id", "imp", "ref", "alt", "dr2", "af"),
      colClasses = c("character", "numeric", "character", "character", "character", "character", "numeric", "numeric"),
      stringsAsFactors = F
    ) 
    n_03 <- nrow(subset(imp, dr2 >= 0.3))
    n_09 <- nrow(subset(imp, dr2 >= 0.9))
    n_03_common <- nrow(subset(imp, dr2 >= 0.3 & af >= 0.01))
    n_09_common <- nrow(subset(imp, dr2 >= 0.9 & af >= 0.01))
    n_common <- nrow(subset(imp, af >= 0.01))
    n_all <- nrow(imp)

     write(
            x = "| Allele frequency filter | DR2 filter | Number of variants |\n|:----|:----|:----|",
            file = md,
            append = T
        )
        write(
            x = paste0("| All variants | None | ", n_all, " |"),
            file = md,
            append = T
        )
        write(
            x = paste0("| All variants | DR2 >= 0.3 | ", n_03, " |"),
            file = md,
            append = T
        )
        write(
            x = paste0("| All variants | DR2 >= 0.9 | ", n_09, " |"),
            file = md,
            append = T
        )
        write(
            x = paste0("| Common variants (AF >= 0.01) | None | ", n_common, " |"),
            file = md,
            append = T
        )

        write(
            x = paste0("| Common variants (AF >= 0.01) | DR2 >= 0.3 | ", n_03_common, " |"),
            file = md,
            append = T
        )
        write(
            x = paste0("| Common variants (AF >= 0.01) | DR2 >= 0.9 | ", n_09_common, " |"),
            file = md,
            append = T
        )

write(
    x = "\n",
    file = md,
    append = T
)
    

    #batr plot of dr2 values
    ggplot(imp, aes(x=af, y=dr2)) + 
        geom_point(alpha=0.3) +
        xlab("Allele frequency") +
        ylab("Dosage R^2") +
        ggtitle(paste0("Chromosome ", chr)) + geom_xsidedensity(
      data = imp,
        mapping = aes(
            x = af,
            y = after_stat(density)
        ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = imp,
      mapping = aes(
        x = after_stat(density),
        y = dr2
      ),
      alpha = 0.8
    ) 
    ggsave(filename = paste0(plot_folder, "imputation_dr2_chr", chr, ".png"), dpi=200, width = 5, height = 5)
    write(
    x = paste0("![](", relative_plot_folder, "imputation_dr2_chr", chr, ".png)"),
    file = md,
    append = T
    )
    af <- ref %>% select(chr, pos, af_ref = af) %>% inner_join(imp %>% select(chr, pos, af_imp = af), by = c("chr", "pos"))


    ggplot(af, aes(x=af_ref, y=af_imp)) + 
        geom_point(alpha=0.3) +
        geom_abline(slope=1, intercept=0, color="red") +
        xlab("Allele frequency in reference panel") +
        ylab("Allele frequency in imputed data") +
        ggtitle(paste0("Chromosome ", chr)) + geom_xsidedensity(
      data = af,
      mapping = aes(
        x = af_ref,
        y = after_stat(density)
      ),
      alpha = 0.8
    ) +
    geom_ysidedensity(
      data = af,
      mapping = aes(
        x = after_stat(density),
        y = af_imp
      ),
      alpha = 0.8
    ) 
    ggsave(filename = paste0(plot_folder, "imputation_af_chr", chr, ".png"), dpi=200, width = 5, height = 5)
    write(
    x = paste0("![](", relative_plot_folder, "imputation_af_chr", chr, ".png)"),
    file = md,
    append = T
    )
}