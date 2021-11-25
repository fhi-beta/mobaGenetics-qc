# ---- 0. Load dependencies

message('Loading script dependencies...\n')


# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(ggplot2, quietly=TRUE)
    library(gridExtra, quietly=TRUE)
    library(grid, quietly=TRUE)
    library(pryr, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

message("printing arguments from snakemake...\n")
print(args)
message("\n")

# improve readability by unpacking args to variables:
input.methylsets = args[1]
input.rgset = args[2]

params.plot_titles = args[3]

output.plots = args[4]

plotTitles <- c('Unprocessed', params.plot_titles)

# ---- 1. Read in methylSet
message("Reading in data from:\n\t", 
        paste0(input.rgset, '\n\t', paste0(input.methylsets, collapse='\n\t'), '\n'))

rgSet <- readRDS(input.rgset)
set.seed(10)
cpg_subsamples <- sample(rownames(rgSet), size=50000)
betaRGSet <- na.omit(getBeta(rgSet[cpg_subsamples, ]))
rm(rgSet)
invisible(gc())

sampleNames <- colnames(betaRGSet)
cpg_subsamples <- rownames(betaRGSet)

methylSets <- lapply(input.methylsets, readRDS)
.getBetaNaOmit <- function(methylSet) na.omit(getBeta(methylSet))[cpg_subsamples, sampleNames]
betaMethylSets <- lapply(methylSets, FUN=.getBetaNaOmit)
rm(methylSets)
invisible(gc())


betaList <- c(list(betaRGSet), betaMethylSets)

rm(betaMethylSets); rm(betaRGSet); 
invisible(gc())


# ----- 2. Rendering Plots and Writing to PDF 
message(paste0("Writing plots to:\n\t", output.plots, '\n'))

pickBetaForSample <- function(betaList, sampleName, plotTitles) {
  n_titles <- length(plotTitles)
  final_df <- data.frame(matrix(ncol=n_titles,nrow=0))
  colnames <- c('beta','method')
  colnames(final_df) <- colnames
  for (i in 1:n_titles) {
    beta <- betaList[[i]][,sampleName]
    method <- rep(plotTitles[i], length(beta))
    df <- cbind(beta, method)
    colnames(df) <- colnames
    final_df <- rbind(final_df, data.frame(beta, method))
  }
  return(final_df)
}


.densityPlotSample <- function(sampleName, betaList, plotTitles) {
  # message(paste0('Plotting: ', sampleName))
  betaDF <- pickBetaForSample(betaList, sampleName, plotTitles)
  density <- ggplot(data=betaDF, aes(x=beta, group=method, color=method)) + geom_density() +
    xlab('Beta Values') + ylab('Density') + ggtitle(sampleName)
  return(density)
}


message('Saving plots to file...\n')
pdf(output.plots, height=6, width=8.5)
#lapply(densities, FUN = function(x) x)
for(i in sampleNames){
print(.densityPlotSample(i, betaList, plotTitles))
}
dev.off()


message("Plotting done!\n\n")
