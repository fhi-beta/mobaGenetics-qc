#!/bin/bash

# Define the base path and other variables
BASE_PATH="/mnt/work/qc_genotypes/pipeOut_dev/2025.01.30/mod6-imputation"
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
BATCHES=("snp001" "snp002" "snp003" "snp007" "snp008" "snp009" "snp010" "snp011" "snp012"
         "snp013" "snp014" "snp015a" "snp015b" "snp016a" "snp016b" "snp017a" "snp017b"
         "snp017c" "snp017d" "snp017e" "snp017f" "snp018a" "snp018b" "snp018c" "snp018de"
         "snp019")

# Function to convert VCF to BCF and index it
convert_to_bcf() {
    local vcf_file=$1
    local bcf_file="$2"
    bcftools view -O b -o "$bcf_file" "$vcf_file" && bcftools index "$bcf_file"
}

export -f convert_to_bcf

# Find all gzipped VCF files and convert them in parallel, save as .bcf files
find "$BASE_PATH" -type f -name "mod6_conform.chr21.conformed.vcf.gz" | \
parallel --jobs $(nproc) convert_to_bcf {} '{=s/.vcf.gz/.bcf/=}'

echo "Conversion completed."