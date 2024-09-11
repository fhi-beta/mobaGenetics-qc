
##
#
# Writes the main documentation.
#
##


# Libraries

library(conflicted)
library(yaml)
library(glue)
library(dplyr)
library(ggplot2)
library(ggside)
library(scico)
library(grid)

# Solve name space conflicts
conflicts_prefer(dplyr::filter)

# General parameters
theme_set(theme_bw(base_size = 14))


# Files

args <- commandArgs(TRUE)

parameters_file <- args[1]

if (!file.exists(parameters_file)) {
  
  stop(glue("Parameters file {parameters_file} not found."))
  
}

output_folder <- args[2]

if (!dir.exists(output_folder)) {

  stop(glue("Output folder {output_folder} not found."))

}

pheno_folder <- args[3]

if (!dir.exists(pheno_folder)) {
  
  stop(glue("Phenotypes folder {pheno_folder} not found."))
  
}

release_folder <- args[4]

# Import config file

config <- read_yaml(parameters_file)


# Phenotypes

get_pheno_file <- function(population) {
  
  if (population == "children") {
    
    return(file.path(pheno_folder, "pheno_child"))

  } else if (population == "mothers") {
    
    return(file.path(pheno_folder, "pheno_mother"))

  } else if (population == "fathers") {
    
    return(file.path(pheno_folder, "pheno_father"))

  } else {
    
    stop(glue("Link to phenotype file not implemented for population {population}."))
    
  }
  
}


# Write readme

readme_file <- file.path(output_folder, "README.md")

write(
  x = glue("# {config$name}\n"),
  file = readme_file,
  append = F
)
write(
  x = glue("{config$description}\n"),
  file = readme_file,
  append = T
)
write(
  x = glue("The documentation corresponds to the analyses version `{config$suffix}`.\n"),
  file = readme_file,
  append = T
)

write(
  x = glue("### Analyses\n"),
  file = readme_file,
  append = T
)

for (analysis_id in names(config$analyses)) {
  
  # Link to analysis-specific md
  
  analysis <- config$analyses[[analysis_id]]
  
  relative_path <- glue("{release_folder}/{analysis_id}.md")
  analysis_md <- file.path(output_folder, relative_path)
  
  write(
    x = glue("- [{analysis$name}]({relative_path}): {analysis$description}\n\n"),
    file = readme_file,
    append = T
  )
  
  # Write analysis-specific md
  
  write(
    x = glue("# {analysis$name}\n"),
    file = analysis_md,
    append = F
  )
  write(
    x = glue("{analysis$description}\n"),
    file = analysis_md,
    append = T
  )
  write(
    x = glue("The variable '{analysis$phenotype}' was extracted from the cleaned phenotypic files (MoBa version '{analysis$pheno_version}', phenotypes QC version '{analysis$pheno_release}'). The following covariates were used: {paste(analysis$covariates, sep = ', ')}.\n\n"),
    file = analysis_md,
    append = T
  )

  # Housekeeping

plots_folder <- file.path(output_folder, release_folder, glue("{analysis_id}_plots"))

if (!dir.exists(plots_folder)) {

  dir.create(plots_folder)

}
  
  # Get phenotype summary statistics
  
  for (population in analysis$populations) {
    
    write(
      x = glue("### {population}\n\n"),
      file = analysis_md,
      append = T
    )
    
    write(
      x = glue("#### Phenotypes\n"),
      file = analysis_md,
      append = T
    )
    
    pheno_values <- read.table(
      file = get_pheno_file(population),
      header = T,
      sep = " "
    )
    
    non_missing_values <- as.numeric(pheno_values[[analysis$phenotype]][!is.na(pheno_values[[analysis$phenotype]])])
    
    n_values <- length(unique(non_missing_values))
    
    if (n_values <= 20) {
      
      write(
        x = glue("| Value | N |\n"),
        file = analysis_md,
        append = T
      )
      write(
        x = glue("| ----- | - |\n"),
        file = analysis_md,
        append = T
      )
      
      for (value in sort(unique(non_missing_values))) {
        
        n <- sum(non_missing_values == value)
        write(
          x = glue("| {value} | {n} |\n"),
          file = analysis_md,
          append = T
        )
        
      }
      
      write(
        x = glue("| Total | {length(non_missing_values)} |\n\n"),
        file = analysis_md,
        append = T
      )
      
      for (covariate in analysis$covariates) {
        
        n_covariates <- length(unique(pheno_values[[covariate]][!is.na(pheno_values[[covariate]])]))
        
        if (n_covariates <= 20) {
          
          plot_data <- data.frame(
            covariate = factor(pheno_values[[covariate]]),
            phenotype = factor(pheno_values[[analysis$phenotype]])
          ) %>% 
            filter(
              !is.na(covariate) & !is.na(phenotype)
            )
          
          plot <- ggplot(
            data = plot_data
          ) +
            geom_bar(
              mapping = aes(
                x = phenotype,
                fill = covariate
              ),
              alpha = 0.8,
              position = "dodge"
            ) +
            scale_x_discrete(
              name = analysis$phenotype_name
            ) +
            scale_fill_manual(
              values = scico(
                n = n_values,
                palette = "batlowK",
                begin = 0.2,
                end = 0.8
              )
            ) +
            facet_grid(
              covariate ~ .
            ) +
            theme(
              legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line = element_line()
            )
          
          png(
            filename = file.path(plots_folder, glue("{population}_{covariate}.png")),
            width = 900,
            height = 600
          )
          grid.draw(plot)
          device <- dev.off()
          
          write(
            x = glue("![]({analysis_id}_plots/{population}_{covariate}.png)"),
            file = analysis_md,
            append = T
          )
          
        } else {
          
          plot_data <- data.frame(
            covariate = pheno_values[[covariate]],
            phenotype = factor(pheno_values[[analysis$phenotype]])
          ) %>% 
            filter(
              !is.na(covariate) & !is.na(phenotype)
            )
          
          plot <- ggplot(
            data = plot_data
          ) +
            geom_violin(
              mapping = aes(
                x = phenotype,
                y = covariate,
                fill = phenotype
              ),
              alpha = 0.5
            ) +
            geom_boxplot(
              mapping = aes(
                x = phenotype,
                y = covariate,
                fill = phenotype
              ),
              width = 0.3,
              outlier.shape = NA,
              alpha = 0.5
            ) +
            geom_xsidehistogram(
              mapping = aes(
                x = phenotype,
                fill = phenotype
              ),
              alpha = 0.5,
              stat = "count"
            ) +
            geom_ysidedensity(
              mapping = aes(
                y = covariate,
                fill = phenotype
              ),
              alpha = 0.5
            ) +
            scale_x_discrete(
              name = analysis$phenotype_name
            ) +
            scale_y_continuous(
              name = covariate
            ) +
            scale_fill_manual(
              values = scico(
                n = n_values,
                palette = "batlowK",
                begin = 0.2,
                end = 0.8
              )
            ) +
            theme(
              ggside.panel.scale = 0.15,
              ggside.panel.grid.minor.y = element_blank(),
              ggside.axis.ticks.x = element_blank(),
              ggside.axis.text.x = element_blank(),
              ggside.panel.background = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line = element_line(),
              ggside.axis.line = element_line()
            )
          
          png(
            filename = file.path(plots_folder, glue("{population}_{covariate}.png")),
            width = 900,
            height = 600
          )
          grid.draw(plot)
          device <- dev.off()
          
          write(
            x = glue("![]({analysis_id}_plots/{population}_{covariate}.png)"),
            file = analysis_md,
            append = T
          )
          
        }
      }
      
    } else {
      
      quantiles <- quantile(
        x = non_missing_values,
        probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
      )
      
      write(
        x = glue("| Quantile | Value |\n"),
        file = analysis_md,
        append = T
      )
      write(
        x = glue("| -------- | ----- |\n"),
        file = analysis_md,
        append = T
      )
      
      for (quantile_i in 1:length(quantiles)) {
        
        write(
          x = glue("| {names(quantiles[quantile_i])} | {unname(quantiles[quantile_i])} |\n"),
          file = analysis_md,
          append = T
        )
        
      }
      
      write(
        x = glue("| N | {length(non_missing_values)} |\n\n"),
        file = analysis_md,
        append = T
      )
      
      for (covariate in analysis$covariates) {
        
        n_covariates <- length(unique(pheno_values[[covariate]][!is.na(pheno_values[[covariate]])]))
        
        if (n_covariates <= 20) {
          
          plot_data <- data.frame(
            covariate = factor(pheno_values[[covariate]]),
            phenotype = pheno_values[[analysis$phenotype]]
          ) %>% 
            filter(
              !is.na(covariate) & !is.na(phenotype)
            )
          
          plot <- ggplot(
            data = plot_data
          ) +
            geom_violin(
              mapping = aes(
                x = covariate,
                y = phenotype,
                fill = covariate
              ),
              alpha = 0.5
            ) +
            geom_boxplot(
              mapping = aes(
                x = covariate,
                y = phenotype,
                fill = covariate
              ),
              width = 0.3,
              outlier.shape = NA,
              alpha = 0.5
            ) +
            geom_xsidehistogram(
              mapping = aes(
                x = covariate,
                fill = covariate
              ),
              alpha = 0.5,
              stat = "count"
            ) +
            geom_ysidedensity(
              mapping = aes(
                y = phenotype,
                fill = covariate
              ),
              alpha = 0.5
            ) +
            scale_x_discrete(
              name = covariate
            ) +
            scale_y_continuous(
              name = analysis$phenotype_name
            ) +
            scale_fill_manual(
              values = scico(
                n = n_covariates,
                palette = "batlowK",
                begin = 0.2,
                end = 0.8
              )
            ) +
            theme(
              ggside.panel.scale = 0.15,
              ggside.panel.grid.minor.y = element_blank(),
              ggside.axis.ticks.x = element_blank(),
              ggside.axis.text.x = element_blank(),
              ggside.panel.background = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line = element_line(),
              ggside.axis.line = element_line()
            )
          
          png(
            filename = file.path(plots_folder, glue("{population}_{covariate}.png")),
            width = 900,
            height = 600
          )
          grid.draw(plot)
          device <- dev.off()
          
          write(
            x = glue("![]({analysis_id}_plots/{population}_{covariate}.png)"),
            file = analysis_md,
            append = T
          )
          
        } else {
          
          plot_data <- data.frame(
            covariate = pheno_values[[covariate]],
            phenotype = pheno_values[[analysis$phenotype]]
          ) %>% 
            filter(
              !is.na(covariate) & !is.na(phenotype)
            )
          
          plot <- ggplot(
            data = plot_data
          ) +
            geom_point(
              mapping = aes(
                x = covariate,
                y = phenotype
              ),
              alpha = 0.2
            ) +
            geom_density2d(
              mapping = aes(
                x = covariate,
                y = phenotype
              ),
              col = "white"
            ) +
            geom_xsidedensity(
              mapping = aes(
                x = phenotype
              ),
              alpha = 0.5,
              fill = "grey80"
            ) +
            geom_ysidedensity(
              mapping = aes(
                y = phenotype
              ),
              alpha = 0.5,
              fill = "grey80"
            ) +
            scale_x_continuous(
              name = covariate
            ) +
            scale_y_continuous(
              name = analysis$phenotype_name
            ) +
            theme(
              ggside.panel.scale = 0.15,
              ggside.panel.grid.minor.y = element_blank(),
              ggside.panel.grid.major.y = element_blank(),
              ggside.axis.ticks = element_blank(),
              ggside.axis.text = element_blank(),
              ggside.panel.background = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line = element_line(),
              ggside.axis.line = element_blank()
            )
          
          png(
            filename = file.path(plots_folder, glue("{population}_{covariate}.png")),
            width = 900,
            height = 600
          )
          grid.draw(plot)
          device <- dev.off()
          
          write(
            x = glue("![]({analysis_id}_plots/{population}_{covariate}.png)"),
            file = analysis_md,
            append = T
          )
          
        }
      }
    }
    
    write(
      x = glue("#### Association results\n"),
      file = analysis_md,
      append = T
    )

    write(
      x = glue("![](regenie/{analysis_id}/figures/pop_{population}_pheno_{analysis$phenotype}_mh.png)\n"),
      file = analysis_md,
      append = T
    )
    
    write(
      x = glue("- [Association results](regenie/{analysis_id}/pop_{population}_pheno_{analysis$phenotype}.md)\n"),
      file = analysis_md,
      append = T
    )
    
    write(
      x = glue("- [Results prior to COJO](regenie_no_cojo/{analysis_id}/pop_{population}_pheno_{analysis$phenotype}.md)\n\n"),
      file = analysis_md,
      append = T
    )
    
  }
}

write(
  x = "#### Errors, questions, and bug report\n",
  file = readme_file,
  append = T
)
write(
  x = "We welcome bug reports, suggestions of improvements, and contributions. Please do not hesitate to open an issue a pull request in the repository.\n\n",
  file = readme_file,
  append = T
)

write(
  x = "#### Code of Conduct\n",
  file = readme_file,
  append = T
)
write(
  x = "As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Convenant](https://www.contributor-covenant.org/) [Code of Conduct for Open Source Projects](CODE_OF_CONDUCT.md).\n",
  file = readme_file,
  append = T
)

write(
  x = "#### License\n",
  file = readme_file,
  append = T
)
write(
  x = "Unless otherwise specified in specific files, this work, the associated source code, and results are released under a [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).\n",
  file = readme_file,
  append = T
)
write(
  x = "![CC BY 4.0 License Logo](https://i.creativecommons.org/l/by/4.0/88x31.png)\n\n",
  file = readme_file,
  append = T
)
