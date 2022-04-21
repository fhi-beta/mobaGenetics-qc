#!/bin/bash

# run 1 local config file. By using this script, you dont have to
# navigate to tmp_config if too large sample size for doing all
# samples in one batch.

cores=12  # number of cores

# define rule to run. if first time, do second_part. if rerun, or the
# list of manual removed samples are known, do third_part
rule=third_part
rule=second_part
rule=first_part
rule=build_rgset
# rule=sample_size_check
# output_path must match whatever is in globalConfig.yaml !
output_path=/mnt/archive/gutorm/methylationRuns
# Sometimes, the default /tmp is tiny. Lets grab something bigger (?)
tmpdir=$output_path/tmp
mkdir -p $tmpdir
export TMPDIR=$tmpdir
# For MoBa - this will typically be where the raw-data are found
# But mostly it is used to set global_config right below
data_path=/mnt/archive/Momics/MomicsPub/Methylation/Datasets
# This global config-file is reused by all sets, in the case multiple
# sets are run
global_config="$data_path"/MetCommon/QC/input/globalConfig.yaml

# only for TSD:
# activate virtual environment enabling snakemake. This can be done
# outside of the script prior to "bash runOne.sh" source
# If done here, remember to deactivate at the end
#.venv/bin/activate

# This directory is harcoded in 0_checkSamplesSize.R so for now we leave it here
tmp_config=tmp_config
rm -rf $tmp_config  # remove if the directory exists from before
mkdir -m775 $tmp_config

if [ $# -eq 0 ]
  then
	echo "Usage: $0 <config-files>"
    	echo "Conlesls tmpfig files are typically .yaml files.."
    	echo "Example: $0 input/met00*.yaml"
        echo "It is assumed that a global config file is set in this script: $global_config"
  else
      echo "$# files were given. Will run snakemake for these config files - and maybe make new ones"
      cp $* $tmp_config
      echo just copied $* to $tmp_config

fi

# define the number of local config files to run snakemake for
len_input=$(ls -1q $tmp_config/ | grep met | grep .yaml | wc -l)
while [[ ${len_input} -gt 0 ]]
do
    set=$(ls -1q $tmp_config/ | grep met | grep .yaml | head -1)   # define config file to run
    echo "Starting configfile $set (of new total $len_input)"
    date
    echo snakemake --core $cores --configfile $global_config $tmp_config/$set --snakefile Snakefile $rule
         snakemake --core $cores --configfile $global_config $tmp_config/$set --snakefile Snakefile $rule
    echo "Finished for $set"
    # The command might have produced more configfiles (if number of sample were to large)
    # (courtesy of 0_checkSampleSize.R (snakemake rule sample_size_check)
    # So we first remove the samplesheet we just finished ...
    mv $tmp_config/$set "$output_path"/${set%".yaml"}  # move configfile to the folder with corresponding results
    # ... and then re-evalute the remaining files to deal with 
    len_input=$(ls -1q $tmp_config/ | grep met | grep .yaml | wc -l)
done

# remove the temporary folder with configs
# rm -rf $tmp_config

# deacivate the virtual environment. Should be uncommented if the activation is done within this script
#deactivate


