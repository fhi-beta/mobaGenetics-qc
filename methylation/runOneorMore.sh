#!/bin/bash

# run 1 local config file. By using this script, you dont have to
# navigate to tmp_config if too large sample size for doing all
# samples in one batch.

cores=12  # number of cores

# define rule to run. if first time, do second_part. if rerun, or the
# list of manual removed samples are known, do third_part
rule=third_part
# output_path must match whatever is in globalConfig.yaml !
output_path=/mnt/archive/gutorm/methylationRuns
# For MoBa - this will typically be where the raw-data are found
# But mostly it is used to set global_config right below
data_path=/mnt/archive/Momics/MomicsPub/Methylation/Datasets
# This global config-file is reused by all sets, in the case multiple
# sets are run
global_config="$data_path"/MetCommon/QC/input/globalConfig.yaml

# activate virtual environment enabling snakemake. This can be done
# outside of the script prior to "bash runOne.sh" source
#.venv/bin/activate

rm -rf tmp_config  # remove if the directory exists from before
mkdir -m775 tmp_config

if [ $# -eq 0 ]
  then
	echo "Usage: $0 <config-files>"
    	echo "Config files are typically .yaml files.."
    	echo "Example: $0 input/met00*.yaml"
        echo "It is assumed that a global config file is set in this script: $global_config"
  else
	echo "$# files where given. Will run snakemake for these config files."
	for file in "$@"
	  do
		cp $file tmp_config
	  done

fi

# define the number of local config files to run snakemake for
len_input=$(ls -1q tmp_config/ | grep met | grep .yaml | wc -l)

while [[ ${len_input} -gt 0 ]]
do
        set=$(ls -1q tmp_config/ | grep met | grep .yaml | head -1)   # define config file to run
        echo "Started set $set"
        date
        echo snakemake --core $cores --configfile $global_config tmp_config/$set --snakefile Snakefile $rule
             snakemake --core $cores --configfile $global_config tmp_config/$set --snakefile Snakefile $rule

	echo "Finished for $set"
        mv tmp_config/$set "$output_path"/${set%".yaml"}  # move configfile to the folder with corresponding results
        len_input=$(ls -1q tmp_config/ | grep met | grep .yaml | wc -l)
done

# remove the temporary folder with configs
rm -rf tmp_config

# deacivate the virtual environment. Should be uncommented if the activation is done within this script
#deactivate


