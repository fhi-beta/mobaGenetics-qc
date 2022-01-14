#!/bin/sh

# move to the output folder, given as input to this script in Snakefile
cd $1

# remove stored intermediate data
rm -r tmp_data
rm -r tmp_results

echo "Done with removing unnecessary data "

