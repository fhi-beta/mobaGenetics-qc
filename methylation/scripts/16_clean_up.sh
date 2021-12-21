#!/bin/sh

# move to the output folder, given as input to this script in Snakefile
cd $1

# remove stored intermediate data
cd processed_data
rm *_1_*
rm *_2_*
rm *_4_*
rm *_5_*
rm *_7_*
rm *_8_*
rm *_10_*
rm *_11_*
rm *filter_probes_without_sex_chr.rds

cd ../qc_results
rm *3_filtered_NA_control_probes.csv
rm *4_detection_pvalues.csv
rm *4_dropped_samples.csv
rm *7_filtered_samples.csv

echo "Done with removing unnecessary data "

