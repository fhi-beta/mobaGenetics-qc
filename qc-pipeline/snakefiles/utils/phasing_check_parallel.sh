#!/bin/bash

# Extract arguments
relations_file=$1
threads=$2
phase_script=$3
bcf_file=$4
output=$5

# Split the input file into smaller files
split_lines=$(wc -l < "$relations_file")
lines_per_file=$((split_lines / threads))
split -l "$lines_per_file" -d --additional-suffix=.txt "$relations_file" "${output}_split_"

# Process each split file in parallel
parallel --jobs "$threads" "python $phase_script --trios {} --bcf $bcf_file --output {}.phase_check --optimize" ::: "${output}_split_"*.txt

# Concatenate results and clean up
first_file=$(ls "${output}_split_"*.phase_check | head -n 1)
head -n 1 "$first_file" > "$output"

for file in "${output}_split_"*.phase_check; do
    if [ "$file" != "$first_file" ]; then
        tail -n +2 "$file" >> "$output"
    fi
done

#rm "${output}_split_"*