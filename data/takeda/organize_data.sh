#! /usr/bin/env bash

set -o nounset -o errexit -o pipefail -x

tar_file=($(ls *.tar))

if [ ! ${#tar_file[@]} -eq 1 ]
then
    echo 'More than one tar file present'
    exit 1
fi

tar_file="${tar_file[0]}"

tar -xf "$tar_file"

files=($(for file in GSM*.gz; do nm="${file/_*/}"; echo "$nm"; done | sort -u))

for file in ${files[@]}
do
    mkdir "$file"
    mv "$file"*.gz "$file"

    ln -s -T $(realpath "$file/$file"*barcodes.tsv.gz) "$file/barcodes.tsv.gz"
    ln -s -T $(realpath "$file/$file"*genes.tsv.gz)    "$file/features.tsv.gz"
    ln -s -T $(realpath "$file/$file"*matrix.mtx.gz)   "$file/matrix.mtx.gz"
done

