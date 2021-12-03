#!/bin/sh

# remove stored intermediate data
cd processed_data
rm *.rds

cd ../qc_results
rm *3_filtered_NA_control_probes.csv
rm *4_detection_pvalues.csv
rm *4_dropped_samples.csv
rm *7_filtered_samples.csv

echo "Done with removing unnecessary data "

