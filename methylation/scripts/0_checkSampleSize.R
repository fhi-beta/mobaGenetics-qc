# print start time of script:
start_time = Sys.time()
message(paste0("The script was started at: \n", start_time, "\n\n"))

# Parse arguments
args = commandArgs(trailingOnly=TRUE) # get character vector of file names, both input, params and output. Must be done like this for logging functionality 

# activate renv if renv = TRUE
if(as.logical(args[1])){
        source("renv/activate.R")
}

args = args[-1]

message("Loading script dependencies...\n")

# Suppres package load messages
suppressMessages({
        library(minfi, quietly=TRUE)
})

message("Printing arguments from snakemake...\n")
print(args)

# improve readability by unpacking args to variables:
input.plateDir = args[1]
params.analysisName = args[2]
params.maxSampleSize = as.numeric(args[3])
output.sampleSheet = args[4]
output.folder = args[5]

message("Read in sample sheet from plate directory... \n")
message("Sample sheet .rds file exist: ", paste0(file.exists(output.sampleSheet)), "\n")

if(file.exists(output.sampleSheet)){
	arrays = readRDS(output.sampleSheet)
}else{
	arrays <- read.metharray.sheet(input.plateDir)
}
message("Dimension of sample sheet: \n")
print(dim(arrays))

# test if the sample composition is of different sampleType (tissues)
if(sum(colnames(arrays) == "SampleType")){
	message("Checking if your samples are a mix of different tissue types.. \n")
        uniq = unique(arrays$SampleType)  # unique items
        if(length(uniq) > 1 & ("B" %in% uniq)){
                stop("Your samples are a mix of children and adults. \n
                       Split your samples into separate data folders and sample sheets. \n")
        }else{
	message("Your samples are not a mix of different tissues. \n")
	}
}

message("Checking sample size compared to max_sample_size... \n")
if(dim(arrays)[1] > params.maxSampleSize){
	# create folder for the current analysis:
	# dir.create(output.folder)  # does not need to be there, as the log file specification in Snakefile creates the folder
	n = dim(arrays)[1]
	start = 1
	nameTag = 1
	# split into several sample sheets and store them in qc_results within their own folder
	while(start < n){
		end = ifelse(start-1+params.maxSampleSize < n, start-1+params.maxSampleSize, n)
		tmp_arrays = arrays[start:end, ]
		if((n-end)/n <= 0.03){
			end = n
			tmp_arrays = arrays[start:end, ]
		}
		if(nameTag < 10){
			nameTag = paste(0, nameTag, sep = "")
		}
		dir.create(paste(output.folder, "-", nameTag, sep = ""))
		tmp_file_name_sheet = paste(output.folder, "-", nameTag, "/tmp_results/", 
					    params.analysisName, "-", nameTag, "_sampleSheet.rds", sep = "")
		# first, create the folders necessary since saveRDS() can not create folders
		dir.create(paste(output.folder, "-", nameTag, "/tmp_results", sep = ""))
		#dir.create(paste(output.folder, "-", nameTag, "/plots", sep = ""))
		#dir.create(paste(output.folder, "-", nameTag, "/processed_data", sep = ""))
		if(file.exists(tmp_file_name_sheet)){
		}else{
			saveRDS(tmp_arrays, file=tmp_file_name_sheet)
		}
		tmp_file_name_config = paste("tmp_config/", params.analysisName, "-", nameTag, ".yaml", sep = "")
		
		file.copy(paste("tmp_config/", params.analysisName, ".yaml", sep = ""), 
			  tmp_file_name_config)
		nameTag = as.numeric(nameTag) + 1
		start = end+1
	}
	stop(paste0("To large sample size compared to the max_sample_size given in global config file: ", params.maxSampleSize, ". The sample sheet will be split into batches by creation of new local config files in tmp_config."))
}else{
	if(file.exists(output.sampleSheet)){
	}else{
	dir.create(output.folder)
	saveRDS(arrays, file = output.sampleSheet)
}
}

end_time = Sys.time()
message(paste0("The script finished at: \n", end_time, "\n"))

message(paste0("The script had a "))
Sys.time() - start_time

message("Checking and saving sample sheet done!!")

