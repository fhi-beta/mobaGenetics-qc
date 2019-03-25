#!/usr/bin/env bash

dir=/media/local-disk2/jjuod/qcrot1/

# Combines sample exclusion lists from one directory
# Args:
# 1. reporting dir/
# 2. current dataset name
function combine_dir {
	for f in $(ls -rt ${1}exclusions*ind.txt)
	do
		echo "working on file $f"
		stage=${f##*exclusions_}
		stage=${stage%_ind.txt}

		awk -v s=$stage -v d=$2 '{print $0, s, d}' ${f} >> ${dir}exclusions_ind_long.txt
	done
	for f in $(ls -rt ${1}exclusions*snp.txt)
	do
		echo "working on file $f"
		stage=${f##*exclusions_}
		stage=${stage%_snp.txt}
		
		awk -v s=$stage -v d=$2 '{print $0, s, d}' ${f} >> ${dir}exclusions_snp_long.txt
	done
}

echo "FID IID FILTER STAGE BATCH" > ${dir}exclusions_ind_long.txt
echo "CHR SNP FILTER STAGE BATCH" > ${dir}exclusions_snp_long.txt

combine_dir ${dir}mod1-data-preparation/rep/ "MOD1"
combine_dir ${dir}both/rep/ "MOD2-BOTH"
combine_dir ${dir}founders/rep/ "MOD2-FOUNDERS"
combine_dir ${dir}offspring/rep/ "MOD2-OFFSPRING"
combine_dir ${dir}mod3-shaping-preparation/rep/ "MOD3"

Rscript lib/format_excl_table.R ${dir} 
