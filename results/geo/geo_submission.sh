#! /bin/usr/env bash

data=('211230_novogene' '220128_A00405_0521_AHH7JVDSX3' '220311_A00405_0541_AHCMMHDSX3' '221028_A00405_0632_AHJTM7DSX5')
dat_dir='../../data'
res_dir=("~/Projects/tamburini-antigen-tracking/results/2022-03-11" "~/Projects/tamburini-antigen-tracking/results/2022-10-28")
sub='20241101_sheridan_scrnaseq'
sums='geo_scrnaseq_md5sums.txt'
meta='geo_scrnaseq_metadata.xlsx'

set -o errexit -o pipefail -o nounset -x

mkdir -p "$sub"

# geo metadata
ln -sr "$meta" "$sub"

# cellranger matrices
tmp=$(mktemp tmp.XXXXX)

nms=(BT_GEX_3-BT_ADT_3
    BT_GEX_4-BT_ADT_4
    BT_GEX_5-BT_ADT_5
    BT_GEX_1-BT_ADT_1
    BT_GEX_2-BT_ADT_2
    BT_GEX_8-BT_ADT_8
    BT_GEX_7-BT_ADT_7
    BT_GEX_6-BT_ADT_6
    stromal_P2-JH300_1
    stromal_P4-JH300_2
    stromal_P4_P2-JH300_3
    DC_P2-JH300_4
    DC_P4-JH300_5
    DC_P4_P2-JH300_6)

for res in ${res_dir[@]}
do
    for nm in ${nms[@]}
    do
        run=$(basename $res)
        dir="$res/$nm"
        mat_dir="$res/$nm/outs"
    
        if [ ! -d "$dir" ]
        then
            continue
        fi

        for file in 'filtered_feature_bc_matrix.h5'
        do
            mat="$sub/${run}_${nm}_$file"
    
            ln -sr "$mat_dir/$file" "$mat"
    
            md5sum "$mat" \
                >> "$tmp"
        done
    done
done

# # Seurat metadata
# for file in *_count_matrix.h5 *_metadata.tsv.gz
# do
#     ln -sr "$file" "$sub"
# 
#     md5sum "$file" \
#         >> "$tmp"
# done

# fastqs
for dat in ${data[@]}
do
    fq=$dat_dir/$dat/*.*q.gz

    ln -sr $fq "$sub"

    cat "$dat_dir/$dat/md5sums.txt" \
        | sort -k2,2 \
        | awk '$2 !~ ".csv$"' \
        >> $tmp
done

# format md5sums
cat $tmp \
    | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
    > $sums

rm $tmp

