#!/bin/bash

# Extract arguments
phase_script=$1
relations_file=$2
bcf_file=$3
output=$4
threads=$5

# Split the input file into smaller files
split_lines=$(wc -l < "$relations_file")
lines_per_file=$((split_lines / threads))
split -l "$lines_per_file" -d --additional-suffix=.txt "$relations_file" "${output}_split_"

# Process each split file in parallel
parallel --jobs "$threads" "python $phase_script --trios {} --bcf $bcf_file --output {}.phase_check --optimize --pf 10000" ::: "${output}_split_"*.txt

# Concatenate results and clean up
first_file=$(ls "${output}_split_"*.phase_check | head -n 1)
cat "$first_file" > "$output"

for file in "${output}_split_"*.phase_check; do
    if [ "$file" != "$first_file" ]; then
        tail -n +2 "$file" >> "$output"
    fi
done

rm "${output}_split_"*