# ---- 0. Load dependencies

# print start time of script:
start_time = Sys.time()
message(paste0("The script was started at: \n", start_time, "\n\n"))

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality

# activate renv if renv = TRUE
if(as.logical(args[1])){
        source("renv/activate.R")
}

args = args[-1]

message('Loading script dependencies...\n')


# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(ggplot2, quietly=TRUE)
    library(gridExtra, quietly=TRUE)
    library(grid, quietly=TRUE)
})


message("Printing arguments from snakemake...\n")
print(args)
message("\n")

# unpack args to improve readability:
input.methylsets = args[1]
input.rgset = args[2]
input.bmiqed_data = args[3]

params.plot_titles = args[4]

output.plots = args[5]

plotTitles <- c('Unprocessed', params.plot_titles, 'BMIQ')

# ---- 1. Read in methylSet
message("Reading in data from:\n\t",
        paste0(input.rgset, '\n\t', paste0(input.methylsets, collapse='\n\t'), '\n\t', input.bmiqed_data))

message("\n\nFiltering the data to random subsample of 50 000 CpGs... \n ")
betaBMIQ <- readRDS(input.bmiqed_data)

#Pick random CpGs  
sampleNames <- colnames(betaBMIQ)
set.seed(10)
cpg_subsample = sample(rownames(betaBMIQ), size = 50000)
betaBMIQ = betaBMIQ[cpg_subsample,]
invisible(gc())

# Read in methylSets and subsample based on the randomly picked CpGs and samples that survived until BMIQ
methylSets <- lapply(input.methylsets, readRDS)
.getBetaNaOmit <- function(methylSet) na.omit(getBeta(methylSet))[cpg_subsample, sampleNames]
betaMethylSets <- lapply(methylSets, FUN=.getBetaNaOmit)
rm(methylSets)
invisible(gc())

# Read in rgSet and subsample based on randomly picked CpGs and samples that survived until BMIQ
rgSet <- readRDS(input.rgset)
betaRGSet <- na.omit(getBeta(subsetByLoci(rgSet, includeLoci = cpg_subsample, keepControls = FALSE, keepSnps = FALSE)))[, sampleNames]
rm(rgSet)
invisible(gc())

# ---- 2. Assemble DataTable of Beta Values for Plotting

message('\n\nBuilding data.table object for plotting...\n')

betaList <- c(list(betaRGSet), betaMethylSets, list(betaBMIQ))

rm(betaBMIQ); rm(betaMethylSets); rm(betaRGSet); 
invisible(gc())


# ----- 3. Rendering Plots and Writing to PDF
# we assume that beta sets and processing methods are in the same order
message(paste0("\n\nWriting plots to:\n\t", output.plots, '\n'))
pickBetaForSample <- function(betaList, sampleName, plotTitles) {
  n_titles <- length(plotTitles)
  final_df <- data.frame(matrix(ncol=n_titles,nrow=0))
  colnames <- c('beta','method')
  colnames(final_df) <- colnames
  # For each of the processing methods pick beta values corresponding to the sampleName and bind them in one DF
  for (i in 1:n_titles) {
    beta <- betaList[[i]][,sampleName]
    method <- rep(plotTitles[i], length(beta))
    df <- cbind(beta, method)
    colnames(df) <- colnames
    final_df <- rbind(final_df, data.frame(beta, method))
  }
  return(final_df)
}


# For given sampleName pick betas for all processing methods and make density plots
.densityPlotSample <- function(sampleName, betaList, plotTitles) {
 # message(paste0('Plotting: ', sampleName))
  betaDF <- pickBetaForSample(betaList, sampleName, plotTitles)
  density <- ggplot(data=betaDF, aes(x=beta, group=method, color=method)) + geom_density() +
    xlab('Beta Values') + ylab('Density') + ggtitle(sampleName)
  return(density)
}


message('\n\nSaving plots to file...\n')
pdf(output.plots, height=6, width=8.5)
for(i in sampleNames){
print(.densityPlotSample(i, betaList, plotTitles))
}
dev.off()

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time


message("Plotting done!\n\n")

