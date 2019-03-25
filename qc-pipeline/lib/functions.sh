# EIGENSTRAT HELPER FUNCTIONS
# 
# Create param file for convertf
# Args:
# 1. Input path + stem (og ped/map fileset)
# 2. Output path + stem
# 3. Output path + filename of convertf param file
function create_convertf_paramfile {
    echo "genotypename:    $1.ped" > $3
    echo "snpname:         $1.pedsnp" >> $3
    echo "indivname:       $1.pedind" >> $3
    echo "outputformat:    EIGENSTRAT" >> $3
    echo "genotypeoutname: $2.geno" >> $3
    echo "snpoutname:      $2.snp" >> $3
    echo "indivoutname:    $2.ind" >> $3
    echo "familynames:     YES" >> $3
}

# Create param file for smartpca
# Args:
# 1. Input path + stem (og ped/map fileset)
# 2. Output path + stem
# 3. Output path + filename of param file
# 4. txt list of anchor populations
function create_smartpca_paramfile {
	echo "genotypename:	$1.geno" > $3
	echo "snpname:		$1.snp" >> $3
	echo "indivname:	$1.ind" >> $3
	echo "evecoutname:	$2.pca.evec" >> $3
	echo "evaloutname:	$2.eval" >> $3
	if [ ! -z ${4:-} ]
	then
		echo "poplistname:	$4" >> $3
	else
		echo "fastmode:		YES" >> $3
	fi
	echo "altnormstyle:	NO" >> $3
	echo "numoutevec:	10" >> $3
	echo "numoutlieriter:	0" >> $3
	echo "numoutlierevec:	10" >> $3
	echo "outliersigmathresh: 6" >> $3
	echo "qtmode:		0" >> $3
}

# Args:
# 1. Input bedset (needs to be pre-pruned)
# 2. Output path + prefix of PCA output
# 3. TMP folder
function pca {
	# Make sure output folder exists
	mkdir -p `dirname $2`

	# Recode
	plink \
	    --bfile $1 \
	    --autosome \
	    --recode \
	    --out $3/pca_recoded_tmp

	# CONVERTF
	# Create convertf compliant fileset (pedsnp and pedind)
	cat $3/pca_recoded_tmp.map > $3/pca_recoded_tmp.pedsnp
	awk '{print $1" "$2" "$3" "$4" "$5" 1"}' $3/pca_recoded_tmp.ped > $3/pca_recoded_tmp.pedind

	create_convertf_paramfile \
		$3/pca_recoded_tmp \
		$3/pca_recoded_eig \
		$3/pca_convertf_params

	PATH=$PATH:bin/eigenstrat/
	convertf -p $3/pca_convertf_params

	# PCA
	create_smartpca_paramfile \
		$3/pca_recoded_eig \
		$2 \
		$3/merge-pca.par

	smartpca -p $3/merge-pca.par
}


# Takes a pre-pruned input dataset, merges with HapMap and runs PCA
# Args:
# 1. Input stem of bedset
# 2. Input stem of hapmap
# 3. Output path + prefix of PCA results
# 4. TMP-folder
function pca_with_hapmap {
	# Make sure output path exists
	mkdir -p `dirname $3`

	# Get list of all markers in input dataset
	awk '{print $2}' $1.bim > $4/pcain_pruned_markerlist

	# Extract markers in input dataset from HapMap data
	$plinklocal \
	    	--bfile $2 \
		--extract $4/pcain_pruned_markerlist \
	    	--make-bed \
	    	--out $4/hapmap_pruned

	# Check to see if multiple chromosomes are seen for single variant
	# Eg. if a variant got new position from b36 to b37
	# Result from this function is a markerset ok for merge and PCA
	awk '{n[$2]++} END{for(s in n) if(n[s]==1) print s}' \
		$4/hapmap_pruned.bim > $4/referencedata.studymarkers.tmp.prune.in
	numsnps_touse=`wc -l $4/referencedata.studymarkers.tmp.prune.in | cut -d ' ' -f1`
	store "NSNP_TOBEUSEDINPCA" ${numsnps_touse}

	# Extract markers from study data
	$plinklocal \
		--bfile $1 \
		--extract $4/referencedata.studymarkers.tmp.prune.in \
		--make-bed \
		--out $4/studydata_refmarkers_ok

	# Extract markers from HapMap data
	$plinklocal \
		--bfile $2 \
		--extract $4/referencedata.studymarkers.tmp.prune.in \
		--make-bed \
		--out $4/refdata_refmarkers_ok

	# merge datasets the first time
	# || true not to trigger error if set -e
	$plinklocal \
		--bfile $4/refdata_refmarkers_ok \
		--allow-no-sex \
		--bmerge $4/studydata_refmarkers_ok \
		--make-bed \
		--out $4/first_merge_pca || true

	# remove markers with mismatching allele codes
	$plinklocal \
		--bfile $4/studydata_refmarkers_ok \
		--exclude $4/first_merge_pca-merge.missnp \
		--make-bed \
		--out $4/studydata

	$plinklocal \
		--bfile $4/refdata_refmarkers_ok \
		--exclude $4/first_merge_pca-merge.missnp \
		--make-bed \
		--out $4/refdata

	#remerge
	$plinklocal \
		--bfile $4/refdata \
		--allow-no-sex \
		--bmerge $4/studydata \
		--make-bed \
		--out $4/merge-pca

	report_nsnps $4/merge-pca "MERGEDFORPCA"

	# recode to ped
	$plinklocal \
		--bfile $4/merge-pca \
		--recode \
		--out $4/merge-pca

	# Create files for convertf
	# Genotype file: use .ped
	# SNP-file:
	cat $4/merge-pca.map > $4/merge-pca.pedsnp
	# Indivfile
	awk -v OFS=' ' '{print $1, $2, $3, $4, $5, $2~/NA/ ? "hapmap" : 1}' $4/merge-pca.ped > $4/merge-pca.pedind
	# Populations
	echo "hapmap" > $4/poplist.txt

	create_convertf_paramfile \
		$4/merge-pca \
		$4/merge-pca-eig \
		$4/merge_pca_convertf_params

	PATH=$PATH:bin/eigenstrat/
	convertf -p $4/merge_pca_convertf_params
	
	create_smartpca_paramfile \
		$4/merge-pca-eig \
		$3 \
		$4/merge-pca.par \
		$4/poplist.txt
	smartpca -p $4/merge-pca.par
	log "Running PCA on merged studydata and HapMap"
}

# # identify family types for report
# # Args:
# # 1. .fam file which to use
# # 2. reporting folder
# function get_family_types {
#         # makes a list with columns N, hasFather, hasMother, numSiblings
#         awk 'FNR==NR{
#                 i[$2]=p[$3]=m[$4]=1; o[$3,$4]++; next
#         }
#         !($2 in p || $2 in m){
#                 o[0,0]=1; print $3 in i, $4 in i, o[$3,$4]
#         }' $1 $1 |
#         sort | uniq -c |
#         awk '{print $1/$4, $2, $3, $4}' > $2/family-structures.txt
#         # make images   
#         Rscript lib/draw-pedigrees-helper.R \
#                 $2/family-structures.txt \
#                 $2/family-structures
# }
# 

# LOGGING HELPERS
# # Clears the log
# function clear_log {
# 	cat /dev/null > ${PIPELINELOG}
# }
# 
# Log message
function log {
	echo $1 >> ${PIPELINELOG}
}

# REPORTING HELPERS
# To clear file for variable storage
function clear_store {
	echo "${STOREDVARS}"
	cat /dev/null > ${STOREDVARS}
	cat /dev/null > ${EXCLUSIONSI}
	cat /dev/null > ${EXCLUSIONSS}
}
# To store simple variable::value combo 
function store {
	echo $1::$2 >> ${STOREDVARS}
}

# Calculates SNP or IND numbers and stores them
# Args:
# 1. binary filestem
# 2. name of current stage
function report_nsnps {
        num=`wc -l $1.bim | cut -f1 -d' '`
	store "NSNP_${2}" $num
}
function report_ninds {
        num=`wc -l $1.fam | cut -f1 -d' '`
	store "NIND_${2}" $num
}

# Logs excluded SNPS or INDS to track them easily
# Args:
# 1. binary filestem BEFORE
# 2. binary filestem AFTER
# 3. name of exclusion filter
function report_lostinds {
	TMP_EXCL=${EXCLUSIONSI}_tmp
	awk -v f=$3 'FNR==NR{a[$1 " " $2]=1; next}
		{if(!($1 " " $2 in a)) print $1, $2, f}' $2.fam $1.fam > $TMP_EXCL
	num=`wc -l $TMP_EXCL | cut -d ' ' -f1`
	cat $TMP_EXCL >> ${EXCLUSIONSI}
	rm $TMP_EXCL
	store "NIND_${3}_LOST" $num
}
function report_lostsnps {
	TMP_EXCL=${EXCLUSIONSS}_tmp
	awk -v f=$3 'FNR==NR{a[$1 " " $2]=1; next}
		{if(!($1 " " $2 in a)) print $1, $2, f}' $2.bim $1.bim > $TMP_EXCL
	num=`wc -l $TMP_EXCL | cut -d ' ' -f1`
	cat $TMP_EXCL >> ${EXCLUSIONSS}
	rm $TMP_EXCL
	store "NSNP_${3}_LOST" $num
}

# Runs PLINK to produce general summary stats
# and R to plot them. To be used after each filter stage.
# Args:
# 1. Filestem of input of the current stage
# 2. Directory for temp storage of summary statistics
# 3. Stage name
# 4. summary statistic (L/I/H/E/M)
# 5. cutoff
# [6. optional: common marker threshold for HET]
function report_plots {
	infile=$1
	stage=$3

	# Plink to generate summaries, and R to plot
	# Final 5 arguments are thresholds for LMISS, IMISS, HET, HWE, MAF.
	case $4 in 
		L) $plinklocal --bfile $infile --missing --out $2/summary-stats
		   ${libdir}/draw-summary-plots.R \
			$2/summary-stats \
			$PLOTDIR \
			$stage \
			$5 NA NA NA NA
		  ;;
		I) $plinklocal --bfile $infile --missing --out $2/summary-stats
		   ${libdir}/draw-summary-plots.R \
			$2/summary-stats \
			$PLOTDIR \
			$stage \
			NA $5 NA NA NA
		  ;;
		H) $plinklocal --bfile $infile --autosome \
			--maf $6 --het --out $2/summary-stats
		   $plinklocal --bfile $infile --autosome \
			--max-maf $6 --het --out $2/summary-stats_rare
		   ${libdir}/draw-summary-plots.R \
			$2/summary-stats \
			$PLOTDIR \
			$stage \
			NA NA $5 NA NA
		  ;;
		E) $plinklocal --bfile $infile --hardy midp --autosome --out $2/summary-stats
		   ${libdir}/draw-summary-plots.R \
			$2/summary-stats \
			$PLOTDIR \
			$stage \
			NA NA  NA $5 NA
		  ;;
		M) echo "MAF summary generation not implemented yet."
		  ;;
		*) echo "----- Flag $4 not recognized, exiting -----"
		   exit 1
		  ;;
	esac
}

# # - Extracts all offspring from the dataset (defined as all samples with
# # non-zero PID or MID).
# # - Selects one random sibling if there are more samples with same PID
# # - PID and MID is set to 0 after exraction.
# # Function returns bedset with extracted offspring only
# # Args:
# # 1. Input dataset
# # 2. Output dataset
# # 3. TMP folder
# function extract_offspring {
# 	# Extracts all offspring (non-founders) from the dataset
# 	tmpfolder=$3
# 
# 	plink \
# 		--bfile $1 \
# 		--filter-nonfounders \
# 		--make-bed \
# 		--out $tmpfolder/all_offspring_only
# 
# 	# Set PID and MID to 0
# 	mv $tmpfolder/all_offspring_only.fam $tmpfolder/all_offspring_only.fam.bak
# 	awk '{print $1" "$2" 0 0 "$5" "$6}' $tmpfolder/all_offspring_only.fam.bak > $tmpfolder/all_offspring_only.fam
# 
# 	# For families with several offspring, randomly extract one offspring from
# 	# each family. The R-script takes a fam file with offspring only and returns
# 	# a list of samples with only one offspring per PID. The seed is necessary
# 	# to reproduce results
# 
# 	lib/select_single_offspring_from_pid_random.R \
# 		$tmpfolder/all_offspring_only.fam \
# 		$tmpfolder/random_single_sibling_from_pid_list \
# 		42
# 
# 	# Extract the obtained list from dataset
# 	plink \
# 		--bfile $tmpfolder/all_offspring_only \
# 		--keep $tmpfolder/random_single_sibling_from_pid_list \
# 		--make-bed \
# 		--out $2
# 
# 	#logging
# 	all_offspring=`wc -l $tmpfolder/all_offspring_only.fam | cut -f1 -d' '`
# 	after_sibremoval=`wc -l $2.fam | cut -f1 -d' '`
# 	log "Extract offspring: ${all_offspring} offspring total. ${after_sibremoval} after removing siblings."
# }

# Manual selection of samples to be extracted as the core sample.
# 1. Eigenstrat PCA output
# 2. HapMap population file
# 3. Output name of filtered samples list
# 4. Output path of saved files
# 5. Std. dev. threshold for inclusion
# 6. Lower limit in first dimension, i.e. left vertical line in the PCA, ignored if NA
# 7. Upper limit in first dimension, i.e. right vertical line in the PCA, ignored if NA
# 8. Lower limit in second dimension, i.e. lower horizontal line in the PCA, ignored if NA
# 9. Upper limit in second dimension, i.e. upper horizontal line in the PCA, ignored if NA
function manual_core_selection {
	Rscript ${libdir}/manual_offspring_selection_moba.R \
		$1 \
		$2 \
		$3 \
		$4 \
		$5 \
		$6 \
		$7 \
		$8 \
		$9

		# Log no of samples extracted
		core_selection=`wc -l $3 | cut -f1 -d' '`
		store "SELECTED_TOKEEP_PCA" ${core_selection} 
		log "Samples in manual extraction list: ${core_selection}"
}


# Removes samples with IBD PIHAT above threshold
# Args:
# 1. Input bedset (function expects pre-pruned input)
# 2. PIHAT threshold
# 3. Output bedset with samples removed
# 4. .genome filestem
function remove_related_above_pihat_threshold {
	# Get samples with PIHAT above 10%
	awk -v pihat="$2" '$10 > pihat && NR>1{print $1, $2, $4, $10}' \
		$4.genome > ${4}_pihat_exclude

	# Remove samples with too high PIHAT
	plink_remove_samples \
		$1 \
		${4}_pihat_exclude \
		$3 \
		"PIHAT_THR"

	#logging
	samples_in=`wc -l $1.fam | cut -f1 -d ' '`
	samples_out=`wc -l $3.fam | cut -f1 -d ' '`
	samples_excluded=`wc -l ${4}_pihat_exclude | cut -f1 -d' '`
	log "Remove by IBD: Samples in: ${samples_in}. Samples out: ${samples_out}. Samples excluded on pihat above ${2}: ${samples_excluded}"
}

# Remove strand ambiguous markers (A/T and C/G SNPs)
# Args:
# 1. input bedset
# 2. output bedset
# 3. TMP-folder
function remove_strand_ambiguous_markers {
	# Get strand ambiguous SNP names
	awk '!($5=="G" && $6=="C") && 
		!($5=="C" && $6=="G") && 
		!($5=="A" && $6=="T") && 
		!($5=="T" && $6=="A"){print $2}' \
		$1.bim > $3/superclean_no_ag_ct_snps

	# Remove strand ambiguous SNPs
	$plinklocal \
	--bfile $1 \
	--extract $3/superclean_no_ag_ct_snps \
	--make-bed \
	--out $2

	#logging
	marker_in=`wc -l $1.bim | cut -f1 -d' '`
	marker_out=`wc -l $2.bim | cut -f1 -d' '`
	ambiguous=`wc -l $3/superclean_no_ag_ct_snps | cut -f1 -d' '`
	report_nsnps $2 "STRAND_nonambiguous"	
	log "Remove ambiguous markers: ${ambiguous} ambiguous. In: ${marker_in}. Out: ${marker_out}."
}

# Function to find the markers present in both of the provided bim files
# Args:
# 1: Input dataset stem 1
# 2: Input dataset stem 2
# 3: Output list of markers in both sets
function get_markers_in_both_bims {
    comm -12 <(cut -f2 $1 | sort) <(cut -f2 $2 | sort) > $3
    present=`wc -l $3 | cut -f1 -d' '`
    log "Markers present in both ${1} and ${2}: ${present}"
}


# # Function compares two input files and output the difference
# # Args
# # 1. file 1
# # 2. file 2
# # 3. output file for diff
# function log_file_diff {
#     comm -3 <(sort $1) <(sort $2) > $3
# }
# 
# 
# # Remove all markers showing discordance where both samples in the duplicate
# # have calls (ie. disregarding one-nc).
# # Args:
# # 1. file with discordance counts per marker (mismatch where both are called)
# # 2. Threshold (>0 get markers with any discordance)
# # 3. Output file containing the discordant markers
# function get_discordant_markers_bothcall {
# 	awk -v threshold="$2" '$2 > threshold { print $1 }' $1 > $3
# 
# 	# logging
# 	num_discordant=`wc -l $3 | cut -f1 -d' '`
# 	log "Duplicate discordance: ${num_discordant} discordant markers with >${2} mismatches detected (both called)"
# 	store "NSNP_DISCORDANT_CALLED" ${num_discordant}
# }
# 
# 
# # Get markers with any mismatch between duplicate pairs. This includes both
# # miscalls and one-nc
# # Args
# # 1. file with discordance counts per marker (any mismatches)
# # 2. Threshold (0 equates to: >0 --> get markers with any discordance)
# # 3. Output file containing the discordant markers
# function get_discordant_markers_mismatch {
# 	awk -v threshold="$2" '$2 > threshold { print $1 }' $1 > $3
# 
# 	# logging
# 	num_discordant=`wc -l $3 | cut -f1 -d' '`
# 	log "Duplicate discordance: ${num_discordant} discordant markers with >${2} mismatches detected (any mismatch)"
# 	store "NSNP_DISCORDANT_ANY" ${num_discordant}
# }
# 
# 

function combine_dir {
     for f in $(ls -rt ${1}exclusions*ind.txt)
     do
         echo "working on file $f"
         stage=${f##*exclusions_}
         stage=${stage%_ind.txt}

         awk -v s=$stage -v d=$2 '{print $0, s, d}' ${f} >> ${3}exclusions_ind_long.txt
     done

     for f in $(ls -rt ${1}exclusions*snp.txt)
     do
         echo "working on file $f"
         stage=${f##*exclusions_}
         stage=${stage%_snp.txt}

         awk -v s=$stage -v d=$2 '{print $0, s, d}' ${f} >> ${3}exclusions_snp_long.txt
     done
}

