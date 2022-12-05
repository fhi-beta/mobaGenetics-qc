#!/usr/bin/env bash
# for debugging, you might want to turn on echo - uncomment next line
set -x
# # Covert TPED to PLINK binary
# # Args:
# # 1. input TPED file
# # 2. filestem for the output
# function plink_tped_to_binary {
#         plink \
#         --tfile $1 \
#         --make-bed \
#         --out $2
# 
#         # logging
#         raw_samples=`wc -l $2.fam | cut -f1 -d' '`
#         raw_markers=`wc -l $2.bim | cut -f1 -d' '`
#         log "Converting GS report to binary. Samples in report: ${raw_samples}, markers: ${raw_markers} markers"
# }
# 
# 
# # Remove duplicates. Requires a file listing samples to remove
# # with FID in first column and IID in second
# # Args:
# # 1. Input plink bed file
# # 2. List of samples to remove
# # 3. Output bedset
# function plink_remove_duplicates {
# 	plink \
# 	--bfile $1 \
# 	--remove $2 \
# 	--make-bed \
# 	--out $3
# 
# 	# logging
# 	samples_input=`wc -l $1.fam | cut -f1 -d' '`
# 	samples_remlist=`wc -l $2 | cut -f1 -d' '`
# 	samples_output=`wc -l $3.fam | cut -f1 -d' '`
# 	report_lostinds $1 $3 "DUPLICATES"
# 	report_ninds $3 "AFTER_DUPLICATES"
# 	store "NIND_TOREMOVE_DUPLICATED" $samples_remlist
# 
# 	log "Remove duplicates: Samples in: ${samples_input}, samples out: ${samples_output}"
# }
# 
# # Update sex info
# # Args:
# # 1. Input plink bed file
# # 2. Sex info file (3 cols: PID IID SEX (1/2 or M/F))
# # 3. Output plink bed file
# function plink_update_sex {
# 	# Update sex info
# 	plink \
# 	--bfile $1 \
# 	--update-sex $2 \
# 	--make-bed \
# 	--out $3
# 
# 	# logging
# 	missing_sex_input=`awk 'BEGIN{a=0} $5==0{a++} END{print a}' $1.fam`
# 	missing_sex_output=`awk 'BEGIN{a=0} $5==0{a++} END{print a}' $3.fam`
# 	nmales=`awk 'BEGIN{a=0} $5==1{a++} END{print a}' $3.fam`
# 	nfemales=`awk 'BEGIN{a=0} $5==2{a++} END{print a}' $3.fam`
# 	samples_to_update=`wc -l $2 | cut -f1 -d' '`
# 	diff_after_update=`comm -13 <(sort $1.fam) <(sort $3.fam) | wc -l | cut -f1 -d' '`
# 
# 	report_ninds $1 "INITIAL"
# 	report_nsnps $1 "INITIAL"
# 	store "NIND_UPDATESEX_MISSINGIN" ${missing_sex_input}
# 	store "NIND_UPDATESEX_MISSINGOUT" ${missing_sex_output}
# 	store "NIND_UPDATESEX_UPDATED" ${diff_after_update}
# 	store "NMALES_INITIAL" ${nmales}
# 	store "NFEMALES_INITIAL" ${nfemales}
# 	report_lostinds $1 $3 "UPDATESEX"
# 	report_ninds $3 "AFTER_UPDATESEX"
# 
# 	log "Update sex: Missing sex input: ${missing_sex_input}, output: ${missing_sex_output}. Samples in update list/updated: ${samples_to_update}/${diff_after_update}"
# 
# }
# 
# # Update sample IDs
# # Args:
# # 1. Input plink binary
# # 2. Sample ID update file
# # 3. Output plink binary
# function plink_update_sample_ids {
# 	plink \
# 	--bfile $1 \
# 	--update-ids $2 \
# 	--make-bed \
# 	--out $3
# 
# 	# logging
# 	not_id_updated=`grep -F -x -f $1.fam $3.fam | wc -l | cut -d' ' -f1`
# 	id_updated=`grep -F -x -v -f $1.fam $3.fam | wc -l | cut -d' ' -f1`
# 
# 	store "NIND_UPDATESAMPLEIDS_UPDATED" ${id_updated}
# 	store "NIND_UPDATESAMPLEIDS_NOTUPDATED" ${not_id_updated}
# 	report_ninds $3 "AFTER_UPDATESAMPLEIDS"
# 	log "Update IDs: number of samples updated: ${id_updated}"
# }
# 
# # Update sample IDs
# # Args:
# # 1. Input plink binary
# # 2. Sample ID update file
# # 3. Output plink binary
# # 4. do not log family types
# function plink_update_parents {
# 	plink \
# 	--bfile $1 \
# 	--update-parents $2 \
# 	--make-bed \
# 	--out $3
# 
# 	# logging
# 	updated=`grep -F -x -v -f $1.fam $3.fam | wc -l | cut -d' ' -f1`
# 	not_updated=`awk 'FNR==NR{i[$2]=1; n=0; next} !($2 in i){n++} END{print n}' $2 $1.fam`
# 
# 	if [ -z "$4" ]; then
# 		# makes a list with columns N, hasFather, hasMother, numSiblings
# 		get_family_types $3.fam $tmp_reporting
# 	fi
# 
# 	store "NIND_UPDATEPARENTS_UPDATED" ${updated}
# 	store "NIND_UPDATEPARENTS_NOTUPDATED" ${not_updated}
# 	# store "SINGLETONS_INITIAL" ${sibships[0]}
# 	# store "SIBS2_INITIAL" ${sibships[1]}
# 	# store "SIBS3PLUS_INITIAL" ${sibships[2]}
# 	log "Update parents: Parental info updated in ${updated} samples, ${not_updated} unchanged."
# }
# 
# # Update allele codes
# # Args:
# # 1. Input plink binary
# # 2. allele update file
# # 3. Output plink binary
# function plink_update_alleles {
# 	plink \
# 	--bfile $1 \
# 	--update-alleles $2 \
# 	--make-bed \
# 	--out $3
# 
# 	# logging
# 	atcg_updated=`grep -F -x -v -f $1.bim $3.bim | wc -l | cut -d ' ' -f1`
# 	atcg_not_updated=`grep -F -x -f $1.bim $3.bim | wc -l | cut -d ' ' -f1`
# 
# 	store "NSNP_ATCG_UPDATED" ${atcg_updated}
# 	store "NSNP_ATCG_NOTUPDATED" ${atcg_not_updated}
# 	grep -F -x -f $1.bim $3.bim >  $tmp_reporting/snps_atcg_notupdated.txt
# 	report_lostsnps $1 $3 "ATCG"
# 	report_nsnps $3 "AFTER_ATCG"
# 
# 	log "Update to ATCG: ${atcg_updated} markers updated, ${atcg_not_updated} unchanged"
# 
# }
# 
# #A function for updating a binary ped file using one of Will's strand files
# # Args:
# #1. The original bed stem (not including file extension suffix)
# #2. The strand file to apply
# #3. The new stem for output
# #4. TMP-folder
# #Result: A new bed file (etc) using the new stem
# function plink_update_build {
# 	#Unpack the parameters into labelled variables
# 	stem=$1
# 	strand_file=$2
# 	outstem=$3
# 	updatetmp=$4
# 
# 	#Cut the strand file into a series of Plink slices
# 	chr_file=$strand_file.chr
# 	pos_file=$strand_file.pos
# 	flip_file=$strand_file.flip
# 	cat $strand_file | cut -f 1,2 > $chr_file
# 	cat $strand_file | cut -f 1,3 > $pos_file
# 	cat $strand_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file
# 
# 	#Because Plink only allows you to update one attribute at a time, we need lots of temp
# 	#Plink files
# 	temp_prefix=TEMP_FILE_XX72262628_
# 	temp1=$updatetmp/$temp_prefix"1"
# 	temp2=$updatetmp/$temp_prefix"2"
# 	temp3=$updatetmp/$temp_prefix"3"
# 
# 	#1. Apply the chr
# 	plink \
# 		--allow-no-sex \
# 		--bfile $stem \
# 		--update-map $chr_file \
# 		--update-chr \
# 		--make-bed \
# 		--out $temp1
# 
# 	#2. Apply the pos
# 	plink \
# 		--allow-no-sex \
# 		--bfile $temp1 \
# 		--update-map $pos_file \
# 		--make-bed \
# 		--out $temp2
# 
# 	#3. Apply the flip
# 	plink \
# 		--noweb \
# 		--allow-no-sex \
# 		--bfile $temp2 \
# 		--flip $flip_file \
# 		--make-bed \
# 		--out $temp3
# 
# 	#4. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file
# 	plink \
# 		--allow-no-sex \
# 		--bfile $temp3 \
# 		--extract $pos_file \
# 		--make-bed \
# 		--out $outstem
# 
# 	#Now delete any temporary artefacts produced
# 	rm -f $updatetmp/$temp_prefix*
# 
# 	# logging
# 	cat "/dev/null" > ${EXCLUSIONS}_snp_changedchr.txt
# 	build_diff=($(awk -v f=${EXCLUSIONS}_snp.txt \
# 			-v fc=${EXCLUSIONS}_snp_changedchr.txt 'BEGIN{bc=bn=0}
# 			FNR==NR{c[$2]=$1; next}
# 			c[$2]!=$1{ bc++ ; print c[$2], $2, "BUILD_CHANGEDCHR" >> fc}
# 			!($2 in c){ bn++ }
# 			{delete c[$2]}
# 			END{for(i in c) print i, c[i], "BUILD" >> f;
# 			print NR-2*FNR, bc, bn}' $stem.bim $outstem.bim))
# 	store "LINES_DIFF"  ${build_diff[0]}
# 	store "CHR_DIFF"  ${build_diff[1]}
# 	store "NAMES_DIFF"  ${build_diff[2]}
# 	report_nsnps $outstem "AFTER_BUILD"
# 	log "Update build: build updated to build 37"
# }
# 
# 
# Exclude markers or regions from list
# Args:
#1. Input bedstem
#2. List of markers OR regions (PLINK format) to exclude
#3. Output bedstem
#4. Exclusion filter
#5. "Do not report" flag
function plink_exclude_markers {
	$plinklocal \
		--bfile $1 \
		--exclude $2 \
		--make-bed \
		--out $3

	# logging
	if [ -z ${5:-} ]; then report_lostsnps $1 $3 "${4}"; fi
	report_nsnps $3 "AFTER_${4}"
}

# Extract markers or region from list
# Args:
# 1. Input bedstem
# 2. List of markers OR regions (PLINK format) to extract
# 3. Output bedstem
# 4. Exclusion filter
# 5. "Do not report" flag
function plink_extract_markers {
	$plinklocal \
		--bfile $1 \
		--extract $2 \
		--make-bed \
		--out $3

	# logging
	marker_in=`wc -l $1.bim | cut -f1 -d' '`
	marker_exlist=`wc -l $2 | cut -f1 -d' '`
	marker_out=`wc -l $3.bim | cut -f1 -d' '`
	log "$4: Markers/regions in extract list ${marker_exlist}. Markers input ${marker_in}, output ${marker_out}"
	
	if [ -z ${5:-} ]; then report_lostsnps $1 $3 "${4}"; fi
	report_nsnps $3 "AFTER_${4}"

}
 
# # Removed chr0 markers from dataset. (Illumina quality control markers.
# # Args:
# # 1. Input bedset
# # 2. Output bedset with chr. 0 removed
# function plink_remove_chr0 {
# 	plink \
# 		--bfile $1 \
# 		--not-chr 0 \
# 		--make-bed \
# 		--out $2
# 
# 	# logging
# 	report_nsnps $2 "AFTER_CHR0"
# 	report_lostsnps $1 $2 "CHR0"
# }
# 
# 
# # Split chromosome X (23) markers in to XY (25) in PAR-region using
# # build 37 (b37) coordinates.
# # Args:
# # 1. Input bedset
# # 2. Output bedset
# function plink_splitx {
# 	plink \
# 		--bfile $1 \
# 		--split-x b37 \
# 		--make-bed \
# 		--out $2
# 
# 	#logging
# 	par_in=`awk 'BEGIN{a=0} $1==25{a++} END{print a}' ${1}.bim`
# 	par_out=`awk 'BEGIN{a=0} $1==25{a++} END{print a}' ${2}.bim`
# 	store "NSNP_PAR_IN" ${par_in}
# 	store "NSNP_PAR_OUT" ${par_out}
# 	log "Split-X: ${par_in} PAR markers in input. ${par_out} markers after split-x"
# }
# 
# # Remove markers below HWE P-value threshold
# function plink_remove_below_hwe_pvalue {
# 	plink \
# 		--bfile $1 \
# 		--hwe $2 \
# 		--make-bed \
# 		--out $3
# 
# 	# logging
# 	markers_input=`wc -l $1.bim | cut -f1 -d' '`
# 	hwe_thr=$2
# 	markers_output=`wc -l $3.bim | cut -f1 -d' '`
# 	log "Remove markers below --hwe p-value ${hwe_thr}. Markers in: ${markers_input}. Markers out: ${markers_output}"
# 
# 	report_lostsnps $1 $3 "HWE${2}"
# }
# 
# Prunes down the markers in the dataset using either
# --indep-pairwise or --indep-pairphase
# Args:
# 1. input dataset
# 2. "indep-pairwise" or "indep-pairphase"
# 3. window-size (variant count)
# 4. step-size (variant count)
# 5. R2-threshold
# 6. output path + bed stem
# 7. TMP-path
function plink_prune_markers {
	$plinklocal \
		--bfile $1 \
		--make-founders \
		--$2 $3 $4 $5 \
		--out $7/prune_markers_tmp

	$plinklocal \
		--bfile $1 \
		--extract $7/prune_markers_tmp.prune.in \
		--make-bed \
		--out $6

	#logging
	markers_in=`wc -l $1.bim | cut -f1 -d' '`
	markers_out=`wc -l $6.bim | cut -f1 -d' '`
	method=$2
	window=$3
	step=$4
	r2=$5
	report_nsnps $6 "AFTER_PRUNE"
	log "Pruning markers by --$2 $3 $4 $5. In: ${markers_in}. Out: ${markers_out}."
}

# Remove markers exceeding X missing genotypes
# Args:
# 1. Input bedset
# 2. geno threshold
# 3. Output bedset
function plink_remove_above_geno {
	$plinklocal \
		--bfile $1 \
		--geno $2 \
		--make-bed \
		--out $3

	# logging
	markers_input=`wc -l $1.bim | cut -f1 -d' '`
	geno_thr=$2
	markers_output=`wc -l $3.bim | cut -f1 -d' '`
	log "Remove markers with --geno above ${geno_thr}. Markers in: ${markers_input}. Markers out: ${markers_output}"
	report_lostsnps $1 $3 "GENO${2}"
	report_nsnps $3 "AFTER_GENO${2}"
}

# # Remove AUTOSOMAL markers exceeding X missing genotypes
# # Args:
# # 1. Input bedset
# # 2. geno threshold
# # 3. Output bedset
# function plink_remove_above_geno_auto {
# 	plink \
# 		--bfile $1 \
# 		--autosome \
# 		--geno $2 \
# 		--make-bed \
# 		--out $3
# 
# 	# logging
# 	markers_input=`wc -l $1.bim | cut -f1 -d' '`
# 	geno_thr=$2
# 	markers_output=`wc -l $3.bim | cut -f1 -d' '`
# 	log "Remove autosomal markers with --geno above ${geno_thr}. Markers in: ${markers_input}. Markers out: ${markers_output}"
# 	report_lostsnps $1 $3 "GENO_AUTO${2}"
# 	report_nsnps $3 "AFTER_GENO_AUTO${2}"
# }

# Removes samples with more than X percent missing genotypes
# Args:
# 1. Input bedset
# 2. mind threshold
# 3. Output bedset
function plink_remove_above_mind {
	$plinklocal \
		--bfile $1 \
		--mind $2 \
		--make-bed \
		--out $3

	# logging
	samples_input=`wc -l $1.fam | cut -f1 -d' '`
	mind_thr=$2
	samples_output=`wc -l $3.fam | cut -f1 -d' '`
	log "Remove samples with --mind above $mind_thr. Samples in: ${samples_input}. Samples out: ${samples_output}"
	report_lostinds $1 $3 "MIND${2}"
	report_ninds $3 "AFTER_MIND${2}"
}

# Keep samples from input list
# Args:
# 1. Input bedset stem
# 2. List of samples to keep (2 cols: PID IID)
# 3. Output bedset stem
# 4. Exclusion filter
function plink_keep_samples {
	$plinklocal \
		--bfile $1 \
		--keep $2 \
		--make-bed \
		--out $3

	# logging
	samples_input=`wc -l $1.fam | cut -f1 -d' '`
	samples_list=`wc -l $2 | cut -f1 -d' '`
	samples_output=`wc -l $3.fam | cut -f1 -d' '`
	report_lostinds $1 $3 "${4}"
	report_ninds $3 "AFTER_${4}"
	store "NIND_TOKEEP_${4}" $samples_list
	log "Samples input: ${samples_input}, keep-list: ${samples_list}, output: ${samples_output}"
}

# Remove samples in input list from binary dataset
# Args:
# 1. Input bedset stem
# 2. List of samples to remove (2 cols: PID IID)
# 3. Output bedset stem
# 4. Exclusion filter
function plink_remove_samples {
	$plinklocal \
		--bfile $1 \
		--remove $2 \
		--make-bed \
		--out $3

	# logging
	samples_input=`wc -l $1.fam | cut -f1 -d' '`
	samples_list=`wc -l $2 | cut -f1 -d' '`
	samples_output=`wc -l $3.fam | cut -f1 -d' '`
	store "NIND_TOREMOVE_${4}" $samples_list
	report_lostinds $1 $3 "${4}"
	report_ninds $3 "AFTER_${4}"

	log "Samples input: ${samples_input}, remove-list: ${samples_list}, output: ${samples_output}"
}

# Runs check-sex command. Returns a list of problematic samples
# Expects no previous pruning of sex markers (function take care of pruning)
# Args:
# 1. Input bedset (all or only sex-markers)
# 2. Female max Fstat (PLINK default 0.20)
# 3. Male min Fstat (PLINK default 0.80)
# 4. Full path to list of samples with PROBLEM according to PLINK check-sex
# 5. TMP folder
function plink_checksex {
        # prune sex markers
        $plinklocal \
                --bfile $1 \
                --chr 23 \
                --indep-pairphase 20000 2000 0.5 \
                --out $5/pruned_sex_markers

        $plinklocal \
                --bfile $1 \
                --extract $5/pruned_sex_markers.prune.in \
                --make-bed \
                --out $5/pruned_for_sexcheck

        # run check-sex flagging problematic samples
        $plinklocal \
                --bfile $5/pruned_for_sexcheck \
                --check-sex $2 $3 \
                --out $5/checksex_report
	
	# draw summary plot
	Rscript ${libdir}/draw-sexcheck-plot.R \
		$5/checksex_report.sexcheck \
		$PLOTDIR/ \
		$3 $2
	
        # extract problem samples from report in list (2 col PID IID)
        awk '$5=="PROBLEM" {print $1, $2}' $5/checksex_report.sexcheck > $4

        # logging total number of samples failing sexcheck
        sex_fail=`awk '$5=="PROBLEM"{print}' $5/checksex_report.sexcheck | wc -l | cut -f1 -d' '`

        if [ $sex_fail -eq 0 ]; then
                log "Check-sex: No samples failed check-sex"
        else
                log "Check-sex: Samples failing check-sex: ${sex_fail}"
        fi
}

# Remove samples with excess heterozygosity.
# Samples outside given std.dev. threshold are removed.
# Only autosomal common markers (passing MAF threshold) are evaluated
# Args:
# 1. Input bedset
# 2. MAF-value for common variants
# 3. SD-threshold for excess het removal
# 4. Output bedset
# 5. TMP folder
function plink_remove_excess_het_common {

        # Common autosomal markers
        $plinklocal \
                --bfile $1 \
                --autosome \
                --maf $2 \
                --het \
                --out $5/common_only_het_tmp

        Rscript ${libdir}/het_fail.R \
                $5/common_only_het_tmp.het \
                $3 \
                $5/common_only_het_fail

	
        # wipe the first row
        tail -n +2 $5/common_only_het_fail > $5/het_fail_common_total

        # Remove samples failing het excess checks
        $plinklocal \
                --bfile $1 \
                --remove $5/het_fail_common_total \
                --make-bed \
                --out $4

        #logging
        n_fail_common_only=`wc -l $5/het_fail_common_total | cut -f1 -d' '`
        samples_input=`wc -l $1.fam | cut -f1 -d' '`
        samples_output=`wc -l $4.fam | cut -f1 -d' '`
        log "Het excess removal - common (MAF > $2) : Samples fail common: ${n_fail_common_only}. Samples in/out ${samples_input} / ${samples_output}"
	report_lostinds $1 $4 "HET${3}"
	report_ninds $4 "AFTER_HET${3}"
}

# Remove samples with excess heterozygosity.
# Samples outside given std.dev. threshold are removed.
# Only autosomal rare markers (under MAF threshold) are evaluated
# Args:
# 1. Input bedset
# 2. MAF-value for rare variant threshold
# 3. SD-threshold for excess het removal
# 4. Output bedset
# 5. TMP folder
function plink_remove_excess_het_rare {

    # Rare autosomal markers
    # Dirty hack 28.10.19 - This can fail, and until we have cleanded up the whole QC, we  accespt failure but lo
    if ! $($plinklocal \
                --bfile $1 \
                --autosome \
                --max-maf $2 \
                --het \
                --out $5/rare_only_het_tmp
	   ); then
	log "Aaargh. plink_remove_excess_het_rare plink failed (no autosoma marker?)"
	return 0  #Lets fake this went ok, ignore plots etc
    fi
        Rscript ${libdir}/het_fail.R \
                $5/rare_only_het_tmp.het \
                $3 \
                $5/rare_only_het_fail

        # wipe the first row
        tail -n +2 $5/rare_only_het_fail > $5/het_fail_rare_total
        tail -n +2 $5/rare_only_het_fail > $5/het_fail_rare_total

        # Remove samples failing het excess checks
        $plinklocal \
                --bfile $1 \
                --remove $5/het_fail_rare_total \
                --make-bed \
                --out $4

        #logging
        n_fail_rare_only=`wc -l $5/het_fail_rare_total | cut -f1 -d' '`
        samples_input=`wc -l $1.fam | cut -f1 -d' '`
        samples_output=`wc -l $4.fam | cut -f1 -d' '`

        log "Het excess removal - rare (MAF < $2): Samples fail rare: ${n_fail_rare_only}. Samples in/out ${samples_input}/${samples_output}"

	report_lostinds $1 $4 "HET_RARE${3}"
	report_ninds $4 "AFTER_HET_RARE${3}"
}


# # Remove markers with excess heterozygosity. Markers are split by common and rare
# # by given MAF threshold. Markers outside given std.dev. threshold are removed.
# # Only autosomal markers are evaluated
# # Args:
# # 1. Input bedset
# # 2. MAF-value for splitting common and rare variants
# # 3. SD-threshold for excess het removal
# # 4. Output bedset
# # 5. TMP folder
# function plink_remove_excess_het_split {
# 
# 	# Common markers
# 	plink \
# 		--bfile $1 \
# 		--autosome \
# 		--maf $2 \
# 		--het \
# 		--out $5/common_het_tmp
# 
# 	Rscript lib/het_fail.R \
# 		$5/common_het_tmp.het \
# 		$3 \
# 		$5/het_fail_common
# 
# 	# Rare markers
# 	plink \
# 	    --bfile $1 \
# 	    --autosome \
# 	    --max-maf $2 \
# 	    --het \
# 	    --out $5/rare_het_tmp
# 
# 	Rscript lib/het_fail.R \
# 	    $5/rare_het_tmp.het \
# 	    $3 \
# 	    $5/het_fail_rare
# 
# 	# Merge het fail files
# 	tail -n +2 $5/het_fail_common > $5/het_fail_total
# 	tail -n +2 $5/het_fail_rare >> $5/het_fail_total
# 
# 	echo "---HET-EXCESS REMOVAL---"
# 	cat $5/het_fail_total
# 
# 	# Remove samples failing het excess checks
# 	plink \
# 		--bfile $1 \
# 		--remove $5/het_fail_total \
# 		--make-bed \
# 		--out $4
# 
# 	#logging
# 	no_fail_common=`tail -n +2 $5/het_fail_common | wc -l | cut -f1 -d' '`
# 	no_fail_rare=`tail -n +2 $5/het_fail_rare | wc -l | cut -f1 -d' '`
# 	samples_input=`wc -l $1.fam | cut -f1 -d' '`
# 	samples_output=`wc -l $4.fam | cut -f1 -d' '`
# 
# 	echo "HET FAIL COMMON"
# 	cat $5/het_fail_common
# 
# 	echo "HET FAIL RARE"
# 	cat $5/het_fail_rare
# 
# 	rm $5/het_fail_common
# 	rm $5/het_fail_rare
# 
# 
# 	log "Het excess removal. Samples fail common: ${no_fail_common}, rare: ${no_fail_rare}. Sampels in/out ${samples_input} / ${samples_output}"
# }
# 
# # Function to remove by Minor Allele Frequency
# # Args:
# # 1 = input bedset
# # 2 = MAF-threshold
# # 3 = output path + filestem
# function plink_remove_below_maf {
# 
# 	plink \
# 		--bfile $1 \
# 		--maf $2 \
# 		--make-bed \
# 		--out $3
# 
# 	#logging
# 	markers_in=`wc -l $1.bim | cut -f1 -d' '`
# 	markers_out=`wc -l $3.bim | cut -f1 -d' '`
# 	maf=$2
# 	log "Remove below MAF ${maf}: Markers in ${markers_in}, markers out ${markers_out}"
# }
# 

# Removes all autosomal markers below Hardy-Weinberg P-value threshold
# Uses the midp flag in PLINK for removal
# Args:
# 1. Input bedset
# 2. p-value threshold
# 3. output path + bedset
# 4. TMP path for storing output
function remove_below_hwe_pval_autosomal {

	$plinklocal \
		--bfile $1 \
		--autosome \
		--hardy midp \
		--out $4/hardy_tmp

	Rscript ${libdir}/hwe_fail.R \
		$4/hardy_tmp.hwe \
		$2 \
		$4/hwe_fail_snplist NA

	$plinklocal \
		--bfile $1 \
		--exclude $4/hwe_fail_snplist \
		--make-bed \
		--out $3

	#logging
	hwe_threshold=$2
	marker_fail=`wc -l $4/hwe_fail_snplist | cut -d ' ' -f1`
	marker_in=`wc -l $1.bim | cut -f1 -d' '`
	marker_out=`wc -l $3.bim | cut -f1 -d' '`
	log "Remove by HWE: Threshold: ${hwe_threshold}. Marker in: ${marker_in}. Marker out: ${marker_out}"

	report_lostsnps $1 $3 "HWE_AUTOSOMAL${2}"
	report_nsnps $3 "AFTER_HWE_AUTOSOMAL${2}"
}

# # Function to filter only founders from the input dataset
# # Args
# # 1. Input bedset
# # 2. Output bedset
# # 3. Tmp folder
# function plink_filter_founders {
# 	awk 'FNR==NR{a[$1]; next} $2 in a' $LIST_FOUNDERS ${1}.fam > $3/founders_tmp.fam
# 
# 	plink \
# 		--bfile $1 \
# 		--keep $3/founders_tmp.fam \
# 		--make-bed \
# 		--out $2
# 
# 	# logging
# 	samples_in=`wc -l $1.fam | cut -f1 -d' '`
# 	samples_out=`wc -l $2.fam | cut -f1 -d' '`
# 	report_ninds $2 "REMAIN_FOUNDERS"
# 	log "Select founders: ${samples_out} founders selected from ${samples_in} samples"
# }
# 
# # Function to filter only non-founders from the input dataset
# # NB Function wipes PID and MID info in the output
# # Args
# # 1. Input bedset
# # 2. Output bedset
# # 3. Tmp folder
# function plink_filter_nonfounders_wipe {
# 	awk 'FNR==NR{a[$1]; next} $2 in a' $LIST_OFFSPRING ${1}.fam > $3/offspring_tmp.fam
# 
# 	plink \
# 		--bfile $1 \
# 		--keep $3/offspring_tmp.fam \
# 		--make-founders \
# 		--make-bed \
# 		--out $2
# 
# 	# logging
# 	samples_in=`wc -l $1.fam | cut -f1 -d' '`
# 	samples_out=`wc -l $2.fam | cut -f1 -d' '`
# 	report_ninds $2	"REMAIN_OFFSPRING"
# 	log "Select non-founders: ${samples_out} non-founders selected from ${samples_in} samples. Pedigree info wiped in output."
# }
# 
# 
# # Removes MID and FID from input dataset
# # Args:
# # 1. Input bedset
# # 2. Output bedset
# function plink_make_founders {
# 	plink \
# 		--bfile $1 \
# 		--make-founders \
# 		--make-bed \
# 		--out $2
# 
# 	# logging
# 	log "Make founders"
# }
# 
# # Merges two binary datasets with defined
# # merge mode - see https://www.cog-genomics.org/plink2/data#merge
# # Args:
# # 1. Input bedset
# # 2. Merge with bedset
# # 3. Merge mode
# # 4. Output
# function plink_bmerge {
# 	plink \
# 		--bfile $1 \
# 		--bmerge $2 \
# 		--merge-mode $3 \
# 		--make-bed \
# 		--out $4
# }
