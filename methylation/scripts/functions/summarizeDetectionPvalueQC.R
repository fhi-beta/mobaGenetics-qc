#' Summarize qc metrics based on probe detection p-values
#' 
#' Takes an RGSet, calculates the proportion and number of failed probes per sample 
#'   at the given `pValue` cut-off then finds the number of samples with proportion failed
#'   probes greater than each of 0.5, 0.1, 0.05, 0.01 and 0.0033.
#' 
#' @param rgSet [`RGChannelSet`] Containing the microarray data from an Illumina EPIC methylation array
#' @param pValue [`numeric`] Float specifying the significance cut-off for calling failed probes
#' 
#' @return [`list`] A list of data tables containing (1) the genomic position by sample detection p-values,
#'   (2) the fraction of failed genomic positions per sample and (3) the number of samples with a proportion
#'   of failed probes greater than each of 0.5, 0.1, 0.05, 0.01 and 0.0033.
#' 
#' @import minfi
#' @import IlluminaHumanMethylationEPICmanifest
#' @import matrixStats
#' @import data.table
#' @export
summarizeDetectionPvalueQC <- function(rgSet, pValue=0.01) {

    resultList <- list()

    message("Computing detection p-values...")

    detectionPvalues <- detectionP(rgSet)
    failed <- detectionPvalues > pValue

    resultList$detectionPvals <- data.table(detectionPvalues, keep.rownames='rownames')

    message('Summarizing probes per sample QC...')

    resultList$sampleQC <- data.table(
        'sample'=colnames(detectionPvalues),
        'fraction_failed_positions'=colMeans(failed),
        'sum_failed_positions'=colSums(failed)
    )

    message('Summarizing number of samples with each proportion of failed probes...')

    resultList$probeQC <- data.table(
        'probes_gt_0.5_samples_failed'=sum(rowMeans(failed) > 0.5),
        'probes_gt_0.1_samples_failed'=sum(rowMeans(failed) > 0.1),
        'probes_gt_0.05_samples_failed'=sum(rowMeans(failed) > 0.05),
        'probes_gt_0.01_samples_failed'=sum(rowMeans(failed) > 0.01),
        'probes_gt_0.0033_samples_failed'=sum(rowMeans(failed) > 0.0033)
    )

    return(resultList)
}
