
##
#
# This script runs some basic QC on the association results and produces the associated documentation.
#
##

#' Returns a ggplot object with the MH for the given association data frame.
#' 
#' @param associationDF the association data frame
#' 
#' @return a ggplot object with the MH
getMh <- function(
    associationDF
) {
  
  # Arrange for plotting
  
  plotDF <- associationDF %>%
    arrange(
      log10p
    ) %>%
    mutate(
      chromosomeNumber = as.numeric(ifelse(chrom == 'X', 23, chrom)),
      x = chromosomeStart[chromosomeNumber] + genpos,
      color = chromosomeNumber %% 2
    ) %>%
    arrange(
      chromosomeNumber, log10p, genpos
    )
  
  maxP <- max(plotDF$log10p)
  
  colors <- c(mhColor1, mhColor2)
  
  plotDF <- plotDF %>% 
    mutate(
      color = factor(color, levels = 0:(length(colors)-1))
    ) %>% 
    arrange(
      color
    )
  
  # Chromosome labels
  
  xLabels <- 1:22
  xLabels[xLabels %% 2 == 0 & xLabels > 17] <- ""
  xLabels <- c(xLabels, "X")
  
  # y axis
  
  maxY <- 5 * ceiling(max(plotDF$log10p / 5))
  maxY <- max(maxY, 10)
  
  yBreaks <- c(0, 5, -log10(5e-8))
  yLabels <- c("", 5, round(-log10(5e-8), digits = 1))
  
  lastBreak <- floor(max(maxY / 10))
  
  if (lastBreak > 0) {
    
    newBreaks <- 10*(1:lastBreak)
    
    while(length(newBreaks) > 3) {
      
      newBreaks <- newBreaks[c(T, F)]
      
    }
    
    yBreaks <- c(yBreaks, newBreaks)
    yLabels <- c(yLabels, round(newBreaks, digits = 1))
    
  }
  
  
  # Build plot
  mhPlot <- ggplot() + 
    geom_hline(
      yintercept = -log10(5e-8), 
      col = "green4", 
      linewidth = 0.3
    ) + 
    geom_point(
      data = plotDF,
      aes(
        x = x, 
        y = log10p, 
        col = color
      ), 
      size = 2
    ) +
    scale_y_continuous(
      name = "p-value [-log10]", 
      breaks = yBreaks, 
      labels = yLabels, 
      expand = expansion(
        mult = c(0, 0.05)
      ), 
      limits = c(0, maxY)
    ) + 
    scale_x_continuous(
      name = "Chromosome", 
      breaks = chromosomeMiddle, 
      labels = xLabels, 
      limits = c(0, genomeLength), 
      expand = expansion(
        mult = 0.01
      )
    ) + 
    scale_color_manual(
      values = colors
    ) + 
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3),
      strip.background = element_rect(
        fill = "grey99"
      )
    )
  
  return(mhPlot)
  
}

#' Returns a ggplot object with the QQ for the given association data frame.
#' 
#' @param associationDF the association data frame
#' 
#' @return a ggplot object with the plot
getQQ <- function(
    associationDF
) {
  
  plotDF <- associationDF %>%
    mutate(
      maf_bin = case_when(
        a1freq < 0.01 ~ "< 1 %",
        a1freq < 0.05 ~ "1-5 %",
        a1freq < 0.1 ~ "5-10 %",
        a1freq >= 0.1 ~ "> 10 %"
      )
    ) %>% 
    group_by(
      maf_bin
    ) %>% 
    arrange(
      log10p
    ) %>% 
    mutate(
      observedLogP = log10p,
      expectedLogP = sort(-log10(ppoints(n = n())))
    ) %>% 
    ungroup()
  
  # Color
  
  plotDF <- plotDF %>% 
    mutate(
      maf_bin = factor(maf_bin, levels = c("< 1 %", "1-5 %", "5-10 %", "> 10 %"))
    ) %>% 
    arrange(
      maf_bin
    )
  
  # Axis
  
  maxValue <- max(plotDF$observedLogP, plotDF$expectedLogP, 10)
  
  
  # Make plot
  qqPlot <- ggplot(
    data = plotDF
  ) + 
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted"
    ) + 
    geom_hline(
      yintercept = -log10(5e-8), 
      col = "green4", 
      linewidth = 0.3
    )  +
    geom_point(
      mapping = aes(
        x = expectedLogP,
        y = observedLogP,
        col = maf_bin
      ),
      size = 2
    ) +
    scale_color_manual(
      values = scico(
        n = 4,
        begin = 0.2,
        end = 0.8
      )
    ) +
    scale_x_continuous(
      name = "Expected p-value [-log10]",
      limits = c(0, maxValue),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "p-value [-log10]", 
      limits = c(0, maxValue),
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(qqPlot)
  
}

#' Returns a ggplot object with the beta plotted against the maf.
#' 
#' @param associationDF the association data frame
#' 
#' @return a ggplot object with the plot
getBetaMaf <- function(
    associationDF
) {
  
  plotDF <- associationDF %>%
    filter(
      is.finite(beta) & is.finite(se) & is.finite(a1freq)
    ) %>% 
    arrange(
      a1freq
    )
  
  
  # Make plot
  betaMafPlot <- ggplot() + 
    geom_hline(
      yintercept = 0
    ) +
    geom_segment(
      data = plotDF,
      mapping = aes(
        x = 100 * a1freq,
        xend = 100 * a1freq,
        y = beta - qnorm(0.975) * se,
        yend = beta + qnorm(0.975) * se
      ),
      col = "grey30",
      linewidth = 0.8
    ) +
    geom_point(
      data = plotDF,
      mapping = aes(
        x = 100 * a1freq,
        y = beta
      ),
      col = "grey20",
      size = 0.8
    ) +
    scale_x_continuous(
      name = "Minor Allele Frequency [%]",
      limits = c(0, 50),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "Effect Size Estimate [95% CI]",
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(betaMafPlot)
  
}

#' Returns a ggplot object with the se plotted against the maf.
#' 
#' @param associationDF the association data frame
#' 
#' @return a ggplot object with the plot
getSeMaf <- function(
    associationDF
) {
  
  plotDF <- associationDF %>%
    filter(
      is.finite(beta) & is.finite(se) & is.finite(a1freq)
    ) %>% 
    arrange(
      a1freq
    )
  
  
  # Make plot
  seMafPlot <- ggplot() + 
    geom_hline(
      yintercept = 0
    ) +
    geom_point(
      data = plotDF,
      mapping = aes(
        x = 100 * a1freq,
        y = se
      ),
      col = "grey20",
      size = 0.8
    ) +
    scale_x_continuous(
      name = "Minor Allele Frequency [%]",
      limits = c(0, 50),
      expand = expansion(
        mult = 0.02
      )
    ) +
    scale_y_continuous(
      name = "Standard Error Estimate",
      expand = expansion(
        mult = 0.02
      )
    )
  
  return(seMafPlot)
  
}


# Libraries

library(conflicted)
library(janitor)
library(stringr)
library(glue)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scico)
library(grid)
library(curl)
library(httr)
library(rjson)


# Namespace conflicts

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Annotation functions

source("gwas/utils/annotation_utils.R")


# Command line input

args <- commandArgs(TRUE)

regenie_output_path <- args[1]
release_folder <- args[2]
analysis_name <- args[3]

# Parameters

theme_set(theme_bw(base_size = 24))


## MH

mhColor1 <- "grey20"
mhColor2 <- "grey40"
bpLimit <- 100000


# Chromosome lengths in GRCh37.p13 (hg19) from Ensembl

chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")
chromosomeLength <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,	135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)
genomeLength <- sum(chromosomeLength)
chromosomeStart <- cumsum(chromosomeLength) - chromosomeLength
chromosomeMiddle <- chromosomeStart + chromosomeLength / 2


# Load the data

regenie_output <- read.table(
  file = regenie_output_path,
  header = T,
  sep = " ",
  stringsAsFactors = F
) %>% 
  clean_names() %>% 
filter(
  is.finite(log10p)
)


# Set up doc files

md_base <- substr(
  x = basename(regenie_output_path),
  start = 7,
  stop = nchar(basename(regenie_output_path)) - 11
)

split <- strsplit(substr(x = md_base, start = 5, stop = 1e9), "_pheno_")[[1]]
population <- split[1]
pheno <- split[2]

doc_folder <- file.path(release_folder, "regenie_no_cojo", analysis_name)
md_file <- file.path(doc_folder, paste0(md_base, ".md"))
figures_folder <- file.path(doc_folder, "figures")
ensembl_folder <- file.path(doc_folder, "ensembl")

# Housekeeping

if (!dir.exists(doc_folder)) {

    dir.create(doc_folder)

}

if (!dir.exists(figures_folder)) {

    dir.create(figures_folder)

}


# Write documentation

write(
  x = glue("## {pheno} in {population}"),
  file = md_file,
  append = F
)

write(
  x = glue("Association results by regenie for {pheno} in {population}.\n"),
  file = md_file,
  append = T
)

if (nrow(regenie_output) == 0) {
  
  write(
    x = glue("**Warning:*** The association results contain no SNP with finite p-value, please check the number of cases vs controls."),
    file = md_file,
    append = T
  )
  
  return()
  
}

write(
  x = glue("### Manhattan"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("{md_base}_mh.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("{md_base}_mh.png"))

write(
  x = glue("![]({relative_figure_path})"),
  file = md_file,
  append = T
)

mh_plot <- getMh(
  associationDF = regenie_output
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(mh_plot)
dummy <- dev.off()

write(
  x = glue("### QQ plot"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("{md_base}_qq.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("{md_base}_qq.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

qq_plot <- getQQ(
  associationDF = regenie_output
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(qq_plot)
dummy <- dev.off()

write(
  x = glue("### Beta vs. Allele Frequency"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("{md_base}_beta_af.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("{md_base}_beta_af.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

beta_af_plot <- getBetaMaf(
  associationDF = regenie_output
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(beta_af_plot)
dummy <- dev.off()

write(
  x = glue("### Standard error vs. Allele Frequency\n"),
  file = md_file,
  append = T
)

absolute_figure_path <- file.path(figures_folder, glue("{md_base}_se_af.png"))
relative_figure_path <- file.path(basename(figures_folder), glue("{md_base}_se_af.png"))

write(
  x = glue("![]({relative_figure_path})\n"),
  file = md_file,
  append = T
)

se_af_plot <- getSeMaf(
  associationDF = regenie_output
)

png(
  filename = absolute_figure_path,
  width = 900,
  height = 600
)
grid.draw(se_af_plot)
dummy <- dev.off()

