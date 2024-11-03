#! /bin/usr/env bash

data=('230623_A00405_0707_AHCK7MDSX7')
dat_dir='../../data'
ref_dir='../../ref'
sub='20241101_sheridan_geomx'
sums='metadata/geo_geomx_md5sums.txt'

set -o errexit -o pipefail -o nounset -x

mkdir -p "$sub"

tmp=$(mktemp tmp.XXXXX)

# RUN THIS ON CLUSTER
# # pkc files
# for file in "$ref_dir/"*.pkc 
# do
#     ln -sr "$file" "$sub"
# 
#     md5sum "$file" \
#         >> "$tmp"
# done
# 
# # GeoMx worksheet and object metadata
# for file in "../geomx/geomx_worksheets/Tamburini_061423_20230619T1930_LabWorksheet.txt" "geomx_metadata.tsv.gz"
# do
#     ln -sr "$file" "$sub"
#     
#     md5sum "$file" \
#         >> "$tmp"
# done
# 
# # fastqs
# for dat in ${data[@]}
# do
#     fq=$dat_dir/$dat/*.fastq.gz
# 
#     ln -sr $fq "$sub"
# 
#     cat "$dat_dir/$dat/md5sums.txt" \
#         | sort -k2,2 \
#         | awk '$2 !~ ".csv$"' \
#         >> $tmp
# done
# 
# # format md5sums
# cat $tmp \
#     | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
#     > $sums

# RUN THIS LOCALLY
# dcc files
# * run this locally since dcc files are not on cluster
for file in ../geomx/counts/*.dcc
do
    ln -sr "$file" "$sub"

    md5sum "$file" \
        >> "$tmp"
done

# format md5sums
cat "$tmp" \
    | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
    >> "$sums"

rm "$tmp"

