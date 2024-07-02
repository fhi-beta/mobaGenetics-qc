#!/bin/bash
# Used to compare results of two different QC-runs. Compares two
# bedsets, comparing .fam and .bim files. For bim files,
# the second column (variant identifier) is ignored (id)
#
# Will sort files before doing a cmp
# The two parameters are the name of the two trunks

if [ ! $# -eq 2 ]
then
    echo Usage: $0 bedsetTrunk1 bedsetTrunk2
    echo Compares  bedsetTrunk1 .fam and .bim files with corresponding bedsetTrunk2
    echo Creates a subdirectory in current directory called sort with
    echo sorted \(and for bim files without column 2\) files
    echo 
    exit 1
fi

# set -x  # For debug/echo
mkdir -p sort
file1=$1;
sortedFile1=sort/$(basename -- "$file1");
file2=$2;
echo comparing $file1.fam  and $file2.fam files with comm -3
# fam files are pure sort
sort $file1.fam > $sortedFile1.fam
sort $file2.fam | comm -3 $sortedFile1.fam -
echo comparing $file1.bim  and $file2.bim files with comm -3
#bim files: Ignore markername and sort order strands ($5 and $6)
awk '{if ($5<$6) {print $1,$3,$4,$5,$6} else {print $1,$3,$4,$6,$5}}' $file1.bim | sort > $sortedFile1.bim
awk '{if ($5<$6) {print $1,$3,$4,$5,$6} else {print $1,$3,$4,$6,$5}}' $file2.bim | sort | comm -3 $sortedFile1.bim -
    
