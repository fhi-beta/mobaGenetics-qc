#!/bin/bash
# Used to compare results of two different QC-runs Checks two
# directories, comparing all plink .fam and .bim files. For bim files,
# se second column (variant identifier) is ignored (id)
#
# Will sort files before doing a cmp
# The two parameters are the name of the directories

if [ ! $# -eq 2 ]
then
    echo Usage: $0 dir1 dir2
    echo Compares all .fam and .bim files found in dir1 and dir2
    echo File names must match \(sorry about that\) and dir1 is master-dir
    echo Creates a subdirectory in current directory called sort with
    echo sorted \(and for bim files without column 2\) files
    echo 
    exit 1
fi

mkdir -p sort
for file in $1/*.{bim,fam};
do
    #echo "$file";
    file1=$1/$(basename -- "$file");
    file2=$2/$(basename -- "$file");
    echo comparing $file1 and $file2 with comm -3
    extension="${file2##*.}"
    if [ $extension == "fam" ]
    then
        sort $file1 > sort/$file
        sort $file2 | comm -3 - sort/$file
    else #bim file: Ignore markername and sort strands
        awk '{if ($5<$6) {print $1,$3,$4,$5,$6} else {print $1,$3,$4,$6,$5}}' $file1 | sort > sort/$file
        awk '{if ($5<$6) {print $1,$3,$4,$5,$6} else {print $1,$3,$4,$6,$5}}' $file2 | sort | comm -3 - sort/$file
    fi
    
done
