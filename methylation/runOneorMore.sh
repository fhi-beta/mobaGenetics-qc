#!/bin/bash

# run 1 local config file. By using this script, you dont have to navigate to tmp_config if too large sample size for doing all samples in one batch.

cores=2  # number of cores

# define rule to run. if first time, do second_part. if rerun, or the list of manual removed samples are known, do third_part
rule=first_part

# activate virtual environment enabling snakemake. This can be done outside of the script prior to "bash runOne.sh"
#source .venv/bin/activate

rm -rf tmp_config  # remove if the directory exists from before
mkdir -m775 tmp_config

if [ $# -eq 0 ]
  then
	echo "No arguments supplied."
    	echo "Will run for all met*.yaml files in the input folder."
    	# copy local config files to tmp_config
	echo "$(ls -1q input/ | grep met | grep .yaml | wc -l) config files detected."
    	
	for file in $(ls -1q input/ | grep met | grep .yaml)
    	  do
        	cp input/$file tmp_config
    	  done
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
        echo snakemake --core $cores --configfile input/globalConfig.yaml tmp_config/$set --snakefile Snakefile $rule
             snakemake --core $cores --configfile input/globalConfig.yaml tmp_config/$set --snakefile Snakefile $rule

        # clean up the output files, removing unecessary files.
	# only do this part if rule=third_part

        if [ "$rule" == "third_part" ]
        then

                echo snakemake --core $cores --configfile input/globalConfig.yaml tmp_config/$set --snakefile Snakefile clean_up
                     snakemake --core $cores --configfile input/globalConfig.yaml tmp_config/$set --snakefile Snakefile clean_up
        fi
	
	echo "Finished for $set"

        mv tmp_config/$set Runs/${set%".yaml"}  # move configfile to the folder with corresponding results

        len_input=$(ls -1q tmp_config/ | grep met | grep .yaml | wc -l)
done

# remove the temporary folder with configs
rm -rf tmp_config

# deacivate the virtual environment. Should be uncommented if the activation is done within this script
#deactivate


