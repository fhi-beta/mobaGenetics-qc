#!/bin/bash

# Runs all sets.
# Environments etc must have been set up correctly, including the global configGlabal.yaml fould found on input
# < more pointers do pre-run documentation>

cores=8 # The number of cores you want to use
# The final snakemake rule to build. Frist time this is run for a set, a manual step needs to be done
# See documentation
rule="third_part" 

for set in input/met*.yaml
do
    echo -n "Started set $set with target rule $rule : "
    date
    echo snakemake -c$cores --configfile $set --snakefile Snakefile $rule
         snakemake -c$cores --configfile $set --snakefile Snakefile $rule
    echo -n "Finished $set with target rule $rule : "
    date
done
