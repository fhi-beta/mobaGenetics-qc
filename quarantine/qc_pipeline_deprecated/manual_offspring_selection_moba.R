#!/usr/bin/env Rscript

#---------------------- Description -------------------------------------------#
#
# This script extracts the population of interest from a PCA file from Eigenstrat.
#
#------------------------------------------------------------------------------#


#---------------------- Attributes --------------------------------------------#

# Import arguments from command line
args <- commandArgs(TRUE)

# File containing the PCA data, must be in the Eigenstrat .pca.evec format.
pcaRawFile <- args[1]
# pcaRawFile = "D:\\projects\\ERC\\PCA\\pruning\\moba12_superclean_core_with_hapmap.pca.evec"

# File containing the hapmap population mapping.
populationMappingFile <- args[2]
# populationMappingFile <- "D:\\projects\\ERC\\PCA\\pruning\\relationships_w_pops_121708.txt"

# File where to save the output
outputFile <- args[3]
# outputFile <- "D:\\projects\\ERC\\PCA\\pruning\\manual\\MoBa12.filtered"

# Folder where to save the plots and other QC files
outputFolder <- args[4]
# outputFolder <- "D:\\projects\\ERC\\PCA\\pruning\\moba12"

# Number of standard deviations to use, ignored if NA
nSTD <- as.numeric(args[5])
# nSTD <- 6

# Lower limit in first dimension, i.e. left vertical line in the PCA, ignored if NA
low1 <- as.numeric(args[6])
# low1 <- NA

# Upper limit in first dimension, i.e. right vertical line in the PCA, ignored if NA
up1 <- as.numeric(args[7])
# up1 <- NA

# Lower limit in second dimension, i.e. lower horizontal line in the PCA, ignored if NA
low2 <- as.numeric(args[8])
# low2 <- NA

# Upper limit in second dimension, i.e. upper horizontal line in the PCA, ignored if NA
up2 <- as.numeric(args[9])
# up2 <- NA

# Value to set for the rejected population
rejectedPopulationName <- "Rejected"

# Output picture width
png_plot_width<-12

# Output picture width
png_plot_height<-10


#---------------------- Dependencies ------------------------------------------#


# install.packages("ggplot2", "~/R", repos="http://cran.uib.no/")
# library(ggplot2, lib.loc = "~/R")
library(ggplot2)
library(class)
library(MASS)
library(RColorBrewer)

#---------------------- Methods ---------------------------------------------#

#' Extracts the id, first and second dimensions of the pca file in a data frame.
#' 
#' @param the file to parse
#' 
#' @return the id, first and second dimensions
parse <- function(inputFile) {
  input <- readLines(inputFile)
  pid <- c()
  iid <- c()
  pca1 <- c()
  pca2 <- c()
  excluded <- c()
  for (lineNumber in 2:length(input)) {
    line <- input[lineNumber]
    column <- 1
    lastIndex <- 1
    space <- T
    for (index in 1:nchar(line)) {
      character <- substring(line, index, index)
      if (character == ' ') {
        if (!space) {
          value <- substring(line, lastIndex, index-1)
          if (column == 1) {
            valueSplit <- strsplit(value, ':')
            tempPid <- valueSplit[[1]][1]
            pid[lineNumber-1] <- tempPid
            tempIid <- valueSplit[[1]][2]
            iid[lineNumber-1] <- tempIid
          } else if (column == 2) {
            numericValue <- as.numeric(value)
            pca1[lineNumber-1] <- numericValue
          } else if (column == 3) {
            numericValue <- as.numeric(value)
            pca2[lineNumber-1] <- numericValue
          }
          excluded[lineNumber-1] <- 0
          column <- column+1
        }
        space <- T
      } else {
        if (space) {
          lastIndex <- index
        }
        space <- F
      }
    }
  }
  result <- data.frame(pid, iid, pca1, pca2, excluded)
  return(result)
}

#' Sets the population according to the ID
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' @param populationMappingFile the file containing the population mapping
#' 
#' @return the updated data frame
setPopulation <- function(pcaResults, populationMappingFile) {
  
  # Read Hapmap populations from file
  populationMapping <- read.table(populationMappingFile, header = T)
  populationMapping <- data.frame(iid=populationMapping$IID, population=populationMapping$population)
  
  # Set HapMap populations from file and the rest to "MoBa"
  pcaResults <- merge(pcaResults, populationMapping, by="iid", all.x=T)
  pcaResults$population <- as.character(pcaResults$population)
  pcaResults$population[is.na(pcaResults$population)] <- "MoBa"
  return(pcaResults)
}

#' Returns the center of a two dimensional distribution given in xy coordinates. 
#' Center is calculated using the median on both dimensions.
#' 
#' @param x the x series
#' @param y the y series
#' 
#' @return the coordinates of the center
getCenter <- function(x, y) {
  return(c(median(x), median(y)))
}


#' Returns the squared distance of points to a center.
#' 
#' @param x the x series
#' @param y the y series
#' @param center the coordinates of the center in c(x0, y0)
#' 
#' @return the coordinates of the center
getDistances <- function(x, y, center) {
  return(sqrt((x-center[1])*(x-center[1])+(y-center[2])*(y-center[2])))
}


#' Returns the maximal distance to the center based on number of estimated standard deviations.
#' 
#' @param distances the distance series
#' @param nStd the number of standard deviations allowed
#' 
#' @return the coordinates of the center
getDistancesLimit <- function(distances, nStd) {
  return(nStd * quantile(distances, 0.68))
}


#' Returns the maximal distance to the center based on percentile.
#' 
#' @param distances the distance series
#' @param q the percentile allowed
#' 
#' @return the coordinates of the center
getDistancesPercentile <- function(distances, q) {
  return(quantile(distances, q))
}

#' Returns the center and limit for std thresholding.
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' @param nSTD the number of standard deviations allowed
#' 
#' @return the filtered data frame
getCenterAndLimit <- function(pcaResults, nSTD) {
  
  # Extract MoBa data
  x <- pcaResults$pca1[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName]
  y <- pcaResults$pca2[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName]
  
  # Get the center of the MoBa population
  center <- getCenter(x, y)
  
  # Get the distances ditribution for the MoBa population and overall population
  mobaDistances <- getDistances(x, y, center)
  allDistances <- getDistances(pcaResults$pca1, pcaResults$pca2, center)
  
  # Get the MoBa distance threshold
  limit <- getDistancesLimit(mobaDistances, nSTD)
  
  return(c(center, limit))
}

#' Filters according to the distance to the center. Sets the population value of the MoBa samples "MoBa_Rejected", and excluded value to 1 to the donors filtered out.
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' @param nSTD the number of standard deviations allowed
#' 
#' @return the filtered data frame
filterStd <- function(pcaResults, nSTD) {
  
  # Get the center and limit
  centerAndLimit <- getCenterAndLimit(pcaResults, nSTD)
  center = c(centerAndLimit[1], centerAndLimit[2])
  limit <- centerAndLimit[3]
  
  # Set excluded to 1 when strictly over the limit
  allDistances <- getDistances(pcaResults$pca1, pcaResults$pca2, center)
  pcaResults$excluded <- ifelse(allDistances > limit, 1, pcaResults$excluded)
  
  # Set population to MoBa_Rejected
  pcaResults$population[pcaResults$population=="MoBa" & pcaResults$excluded > 0] <- rejectedPopulationName
  
  return(pcaResults)
}

#' Filters according to an arbitrary value.
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' @param up boolean indicating whether it is an upper or lower limit
#' @param dimension the dimension number, 1 or 2
#' @param value the value to use as thresold
#' 
#' @return the filtered data frame
filterValue <- function(pcaResults, up, dimension, value) {
  
  # Set excluded to 1 when strictly over the value
  if (dimension == 1) {
    if (up) {
      pcaResults$excluded <- ifelse(pcaResults$pca1 > value, 1, pcaResults$excluded)
    } else {
      pcaResults$excluded <- ifelse(pcaResults$pca1 < value, 1, pcaResults$excluded)
    }
  } else {
    if (up) {
      pcaResults$excluded <- ifelse(pcaResults$pca2 > value, 1, pcaResults$excluded)
    } else {
      pcaResults$excluded <- ifelse(pcaResults$pca2 < value, 1, pcaResults$excluded)
    }
  }
  
  # Set population to MoBa_Rejected
  pcaResults$population[pcaResults$population=="MoBa" & pcaResults$excluded > 0] <- rejectedPopulationName
  
  return(pcaResults)
}

#' Infers the population of the moba population based on the HapMap
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' 
#' @return a data frame containing the moba data with inferren populations
inferPopulation <- function(pcaResults) {
  
  # Make a training set from the hapmap
  trainingMatrix <- data.frame(
    pca1=pcaResults$pca1[pcaResults$population != "MoBa" & pcaResults$population != rejectedPopulationName], 
    pca2=pcaResults$pca2[pcaResults$population != "MoBa" & pcaResults$population != rejectedPopulationName]
  )
  hapmapPopulation <- pcaResults$population[pcaResults$population != "MoBa" & pcaResults$population != rejectedPopulationName]
  
  # get Moba data
  testMatrix <- data.frame(
    pca1=pcaResults$pca1[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName], 
    pca2=pcaResults$pca2[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName]
  )
  
  # Infer Moba population
  inferredPopulations <- knn(train=trainingMatrix, test=testMatrix, cl=hapmapPopulation, k=1, l=0)
  
  # add inferred population to the results
  mobaIds <- pcaResults$iid[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName]
  mobaInferred <- data.frame(iid=mobaIds, pca1=testMatrix$pca1, pca2=testMatrix$pca2, inferredPopulation=inferredPopulations)
  
  return(mobaInferred)
}

#' Draws a circle (taken from http://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2)
#' 
#' @param x0 the center abscissa
#' @param y0 the center ordinate
#' @param r the radius
#' @param n the number of points
#' 
#' @return a data frame containing the points making the circle
getCircle <- function(x0, y0, r, n){
  tt <- seq(0,2*pi,length.out = n)
  xx <- x0 + r * cos(tt)
  yy <- y0 + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#' Plots the two first dimensions with the thresholds
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' 
#' @return the plot
drawGeneralPlot <- function(pcaResults) {
  # Color according to the populations
  labels <- sort(unique(pcaResults$population))
  colfunc <- colorRampPalette(c("Red", "Blue"))
  colors <- colfunc(length(labels))
  colors[match("CEU", labels)] <- "Black"
  colors[match("MoBa", labels)] <- "Green"
  colors[match(rejectedPopulationName, labels)] <- "Grey"
  colors[length(colors) + 1] <- "Black"
  
  # get the range
  min <- min(pcaResults$pca1, pcaResults$pca2)
  max <- max(pcaResults$pca1, pcaResults$pca2)
  range <- c(min, max) 
  
  # Get the title
  nExcluded <- length(pcaResults$population[pcaResults$population == rejectedPopulationName])
  total <- length(pcaResults$population[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName])
  kept <- total - nExcluded
  title <- paste("Excluded: ", nExcluded, " (", round(100 * nExcluded/total), " %); Retained: ", kept, sep="")
  
  # Make plot
  plot <- ggplot()
  plot <- plot + geom_point(data=pcaResults, aes(x=pcaResults$pca1, y=pcaResults$pca2, col=pcaResults$population), size=0.7, alpha=0.5)
  
  if (!is.na(nSTD)) {
    
    # Get the center and limit
    centerAndLimit <- getCenterAndLimit(pcaResults, nSTD)
    
    # Draw circle
    circle <- getCircle(centerAndLimit[1], centerAndLimit[2], centerAndLimit[3], 100)
    plot <- plot + geom_path(data=circle, mapping=aes(x=x, y=y))
  }
  if (!is.na(low1)) {
   plot <- plot + geom_vline(xintercept=low1)
  }
  if (!is.na(up1)) {
    plot <- plot + geom_vline(xintercept=up1)
  }
  if (!is.na(low2)) {
    plot <- plot + geom_hline(yintercept=low2)
  }
  if (!is.na(up2)) {
    plot <- plot + geom_hline(yintercept=up2)
  }
  
  plot <- plot + 
    theme(plot.title=element_text(size=12), legend.key.size=grid::unit(5, "mm")) +
    labs(x="PC1", y="PC2") + 
    guides(colour=guide_legend(title="Population", override.aes=list(size=2))) + 
    scale_x_continuous(breaks=pretty(range, n=5), limits=range) +
    scale_y_continuous(breaks=pretty(range, n=5), limits=range) + 
    scale_colour_manual(values = colors)
  
  return(plot)
}

#' Plots the two first dimensions with the thresholds with zoom on the retained population
#' 
#' @param pcaResults the results of the PCA as parsed by the parse function
#' 
#' @return the plot
drawZoomedPlot <- function(pcaResults) {
  
  # Extracted the pruned data
  pca1 <- pcaResults$pca1[pcaResults$excluded == 0]
  pca2 <- pcaResults$pca2[pcaResults$excluded == 0]
  prunedPopulation <- pcaResults$population[pcaResults$excluded == 0]
  pruned <- data.frame(pca1, pca2, prunedPopulation)
  
  # Color according to the populations
  labels <- sort(unique(prunedPopulation))
  colfunc <- colorRampPalette(c("Red", "Blue"))
  colors <- colfunc(length(labels))
  colors[match("CEU", labels)] <- "Black"
  colors[match("MoBa", labels)] <- "Green"
  colors[match(rejectedPopulationName, labels)] <- "Grey"
  colors[length(colors) + 1] <- "Black"
  
  # Get the title
  nExcluded <- length(pcaResults$population[pcaResults$population == rejectedPopulationName])
  total <- length(pcaResults$population[pcaResults$population == "MoBa" | pcaResults$population == rejectedPopulationName])
  kept = total - nExcluded
  title = paste("Excluded: ", nExcluded, " (", round(100 * nExcluded/total), " %); Retained: ", kept, " - Zoom", sep="")
  
  # Make plot
  plot <- ggplot()
  plot <- plot + geom_point(data=pruned, aes(x=pca1, y=pca2, col=prunedPopulation), size=2, alpha=0.5)
  
  if (!is.na(nSTD)) {
    
    # Get the center and limit
    centerAndLimit <- getCenterAndLimit(pcaResults, nSTD)
    
    # Draw circle
    circle <- getCircle(centerAndLimit[1], centerAndLimit[2], centerAndLimit[3], 100)
    plot <- plot + geom_path(data=circle, mapping=aes(x=x, y=y)) 
  }
  if (!is.na(low1)) {
    plot <- plot + geom_vline(xintercept=low1)
  }
  if (!is.na(up1)) {
    plot <- plot + geom_vline(xintercept=up1)
  }
  if (!is.na(low2)) {
    plot <- plot + geom_hline(yintercept=low2)
  }
  if (!is.na(up2)) {
    plot <- plot + geom_hline(yintercept=up2)
  }
  
  plot <- plot + 
    theme(plot.title=element_text(size=12), legend.key.size=grid::unit(5, "mm")) +
    labs(x="PC1", y="PC2") + 
    guides(colour=guide_legend(title="Population")) + 
    scale_x_continuous(breaks=pretty(pruned$pca1, n=5), limits=c(min(pruned$pca1), max(pruned$pca1))) +
    scale_y_continuous(breaks=pretty(pruned$pca2, n=5), limits=c(min(pruned$pca2), max(pruned$pca2))) + 
    scale_colour_manual(values = colors)
  return(plot)
}

#' Draws the population inference plot
#' 
#' @param inferrenceResults the data frame of the population inferrence
#' 
#' @return the plot
drawInferredPlot <- function(inferrenceResults) {
  
  # Color according to the populations
  labels <- sort(unique(inferrenceResults$inferredPopulation))
  colors <- rainbow(length(labels))
  colors[match("CEU", labels)] <- "Black"
  
  # Get the title
  nCEU <- length(inferrenceResults$inferredPopulation[inferrenceResults$inferredPopulation == "CEU"])
  total <- length(inferrenceResults$inferredPopulation)
  title = paste("MoBa Inferred - CEU: ", nCEU, " (", round(100 * nCEU/total), " %)", sep="")
  
  # Make plot
  plot <- ggplot()
  plot <- plot + geom_point(data=inferrenceResults, aes(x=pca1, y=pca2, col=inferredPopulation), size=0.7, alpha=0.5)
  
  plot <- plot + 
    ggtitle(title) + 
    theme(plot.title=element_text(size=12), legend.key.size=grid::unit(5, "mm")) +
    labs(x="PC1", y="PC2") + 
    guides(colour=guide_legend(title="Population", override.aes=list(size=2))) + 
    scale_x_continuous(breaks=pretty(inferrenceResults$pca1, n=5), limits=c(min(inferrenceResults$pca1), max(inferrenceResults$pca1))) +
    scale_y_continuous(breaks=pretty(inferrenceResults$pca2, n=5), limits=c(min(inferrenceResults$pca2), max(inferrenceResults$pca2))) + 
    scale_colour_manual(values = colors)
  
  return(plot)
}

#' Draws the population inference plot
#' 
#' @param inferrenceResults the data frame of the population inferrence
#' 
#' @return the plot
drawInferredPlotZoomed <- function(inferrenceResults, pcaResults) {
  
  # Extracted the pruned data
  pca1 <- pcaResults$pca1[pcaResults$excluded == 0]
  pca2 <- pcaResults$pca2[pcaResults$excluded == 0]
  prunedPopulation <- pcaResults$population[pcaResults$excluded == 0]
  pruned <- data.frame(pca1, pca2, prunedPopulation)
  
  # Color according to the populations
  labels <- sort(unique(inferrenceResults$inferredPopulation))
  colors <- rainbow(length(labels))
  colors[match("CEU", labels)] <- "Black"
  
  # Get the title
  nCEU <- length(inferrenceResults$inferredPopulation[inferrenceResults$inferredPopulation == "CEU"])
  total <- length(inferrenceResults$inferredPopulation)
  title = paste("MoBa Inferred - CEU: ", nCEU, " (", round(100 * nCEU/total), " %) - Zoom", sep="")
  
  # Make plot
  plot <- ggplot()
  plot <- plot + geom_point(data=inferrenceResults, aes(x=pca1, y=pca2, col=inferredPopulation), size=2, alpha=0.5)
  
  range1 <- c(min(pruned$pca1), max(pruned$pca1))
  range2 <- c(min(pruned$pca2), max(pruned$pca2))
  
  plot <- plot + 
    ggtitle(title) + 
    theme(plot.title=element_text(size=12), legend.key.size=grid::unit(5, "mm")) +
    labs(x="PC1", y="PC2") + 
    guides(colour=guide_legend(title="Population")) + 
    scale_x_continuous(breaks=pretty(range1, n=5), limits=range1) +
    scale_y_continuous(breaks=pretty(range2, n=5), limits=range2) + 
    scale_colour_manual(values = colors)
  
  return(plot)
}


#---------------------- Main Method ---------------------------------------------#

# Parse the input files
pcaResults <- parse(pcaRawFile)
pcaResults <- setPopulation(pcaResults, populationMappingFile)

# Extract the coordinates of the MoBa population
x <- pcaResults$pca1[pcaResults$population == "MoBa"]
y <- pcaResults$pca2[pcaResults$population == "MoBa"]

# If needed, filter according on the STD
if (!is.na(nSTD)) {
  pcaResults <- filterStd(pcaResults, nSTD)
}

# If needed, filter the first dimension on an arbitrary upper value
if (!is.na(up1)) {
  pcaResults <- filterValue(pcaResults, T, 1, up1)
}

# If needed, filter the first dimension on an arbitrary lower value
if (!is.na(low1)) {
  pcaResults <- filterValue(pcaResults, F, 1, low1)
}

# If needed, filter the second dimension on an arbitrary upper value
if (!is.na(up2)) {
  pcaResults <- filterValue(pcaResults, T, 2, up2)
}

# If needed, filter the second dimension on an arbitrary lower value
if (!is.na(low2)) {
  pcaResults <- filterValue(pcaResults, F, 2, low2)
}

# Infer the moba populations
mobaPopulations <- inferPopulation(pcaResults)

# Write all values to the output
write.table(pcaResults, file = file.path(outputFolder, "pca.txt"), sep="\t", quote=F, row.names=F)
write.table(mobaPopulations, file = file.path(outputFolder, "moba_populations.txt"), sep="\t", quote=F, row.names=F)

# Export retained ids to the output file
result <- data.frame(pid=pcaResults$pid[pcaResults$excluded == 0 & pcaResults$population == "MoBa"], iid=pcaResults$iid[pcaResults$excluded == 0 & pcaResults$population == "MoBa"])
write.table(result, file = outputFile, sep="\t", quote=F, row.names=F, col.names=F)

# Make the plots
generalPlot <- drawGeneralPlot(pcaResults)
zoomedPlot <- drawZoomedPlot(pcaResults)
inferencePlot <- drawInferredPlot(mobaPopulations)
inferencePlotZoomed <- drawInferredPlotZoomed(mobaPopulations, pcaResults)
densityPlotAll <- kde2d(mobaPopulations$pca1[], mobaPopulations$pca2[], n=200)
densityPlotCeu <- kde2d(mobaPopulations$pca1[mobaPopulations$inferredPopulation=="CEU"], mobaPopulations$pca2[mobaPopulations$inferredPopulation=="CEU"], n=200)
densityPlotNonCeu <- kde2d(mobaPopulations$pca1[mobaPopulations$inferredPopulation!="CEU"], mobaPopulations$pca2[mobaPopulations$inferredPopulation!="CEU"], n=200)

# Export to pdf
ggsave(filename=file.path(outputFolder, "pca_pruning_general.png"),
	width=png_plot_width, height=png_plot_height, units="cm", plot=generalPlot)
ggsave(filename=file.path(outputFolder, "pca_pruning_zoomed.png"), 
	width=png_plot_width, height=png_plot_height, units="cm", plot=zoomedPlot)
ggsave(filename=file.path(outputFolder, "pca_pruning_inference.png"),
 	width=png_plot_width, height=png_plot_height, units="cm", plot=inferencePlot)
ggsave(filename=file.path(outputFolder, "pca_pruning_inferencezoomed.png"),
 	width=png_plot_width, height=png_plot_height, units="cm", plot=inferencePlotZoomed)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
colorRamp <- rf(32)
pdf(file.path(outputFolder, "pca_remaining_images.pdf"))
image(densityPlotAll, col=colorRamp, xlab="PCA1", ylab="PCA2", main="MoBa12")
image(densityPlotCeu, col=colorRamp, xlab="PCA1", ylab="PCA2", main="MoBa12 inferred CEU")
image(densityPlotNonCeu, col=colorRamp, xlab="PCA1", ylab="PCA2", main="MoBa12 inferred non-CEU")
dev.off()
