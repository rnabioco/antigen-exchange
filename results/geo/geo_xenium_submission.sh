#! /bin/usr/env bash

dat_dir='../xenium/20240322__194506__032224_Tamburini_Run_1'
sub='20241101_sheridan_xenium'
sums='metadata/geo_xenium_md5sums.txt'

set -o errexit -o pipefail -o nounset -x

mkdir -p "$sub"

# raw data
tmp=$(mktemp tmp.XXXXX)

dirs="$dat_dir/"*__?__*

for dir in ${dirs[@]}
do
    dir=$(basename "$dir")

    for file in 'transcripts.csv.gz' "$dat_dir/$dir/"*.ome.tif
    do
        file=$(basename "$file")
        nm=$(echo "$dir" | cut -d '_' -f 3-5)
        nm="${nm}_${file}"

        ln -sr "$dat_dir/$dir/$file" "$sub/$nm"

        md5sum "$sub/$nm" \
            >> "$tmp"
    done
done

# format md5sums
cat $tmp \
    | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
    > $sums

rm $tmp

