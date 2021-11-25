# ---- 0. Load dependencies

message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(wateRmelon, quietly=TRUE)
    library(RPMM, quietly=TRUE)
    library(pryr, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 
message("Printing arguments from snakemake...\n")
print(args)
message("\n")
# unpack args to improve readability:
input.ratioset = args[1]
output.bmiqed = args[2]
output.raw = args[3]

# ---- 1. load GRSet
message(paste0('Loading GRSet no sex chromosomes and getting beta values... \n'))


grSet = readRDS(input.ratioset)

data = getBeta(grSet)

# ---- 2. store NA values for later
message(paste0('Dimension of data table, rows: ', dim(data)[1], ' cols: ', dim(data)[2]))
sum(is.na(data))
NApositionMatrix = is.na(data)

for(i in 1:dim(data)[1]){
	if(sum(is.na(data[i,])) > 0){
 		data[i,is.na(data[i,])] = median(data[i,], na.rm = T)}
}
sum(is.na(data))

# ---- 3. load annotation file:
message('\nReading annotation...')
annotation = getAnnotation(grSet)
rm(grSet)
invisible(gc())

# ---- 4. do bmiq
message(paste0('\nDoing BMIQ...'))
cg = intersect(rownames(data), rownames(annotation))

# source the adjusted BMIQseed function ensuring same set of CpGs i fitting for each sample:
source("scripts/functions/BMIQseed.R")

# normalizes beta values based on probe design. Projects probe II distribution onto probe I
beta.autosomal.BMIQ.list <- apply(data[cg,], 2,BMIQseed, design.v = ifelse(annotation[cg, "Type"] == "I", 1, 2), plots = FALSE, nfit=20000)

# store the normalized beta values
beta.autosomal.BMIQ       <- sapply(beta.autosomal.BMIQ.list, function(x)x$nbeta)
# create possibility of storing the parameters for each sample
beta.autosomal.BMIQ.param <- sapply(beta.autosomal.BMIQ.list, function(x)c(x$av1,x$av2, x$th1,x$th2))

message(paste0('\nThe colnames of the data before and after bmiq is identical: ', identical(colnames(data), colnames(beta.autosomal.BMIQ))))

# ---- 5. Make those CpGs that were NA NA again
for(i in 1:dim(beta.autosomal.BMIQ)[2]){
	beta.autosomal.BMIQ[is.na(NApositionMatrix[,i]),i] = NA
}

message(paste0('Saving beta values prior to BMIQ to ', output.raw))
saveRDS(data, file = output.raw)

message(paste0('Saving bmiqed beta values to ', output.raw))
saveRDS(beta.autosomal.BMIQ, file = output.bmiqed)

message('BMIQ done!')

