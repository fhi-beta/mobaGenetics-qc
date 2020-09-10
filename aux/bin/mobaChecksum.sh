#!/bin/sh
if [  $# -eq 0 ]
then
    echo Usage: $0 file-or-directory ...
    echo "Runs a checksum/hash (sha256sum) on contents of all file/directories (recursively)"
    echo Outputs parameters and hash
    echo "Note that path matters here: dir is not the same as ./dir"
    exit 1
fi
hash=sha256sum

echo $hash of $@
echo To recreate this output, use $0
find $@ -type f -print0 | sort -z | xargs -0 $hash | $hash
