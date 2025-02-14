#!/usr/bin/env bash

set -x

in_fn=$1
out_dir=$2
gffread="/ccb/salz4-3/hji20/hg002-q100-annotation/tools/gffread/gffread"

mkdir -p $out_dir
# sort the reference annotation
$gffread -O -F --keep-exon-attrs $in_fn > "${out_dir}/chm13.sorted.gff"
# remove rRNA, VDJ, and additional gene copies & introduce unique IDs to all features
./clean_chm13.py "${out_dir}/chm13.sorted.gff" $out_dir gff gff

in_fn_2="${out_dir}/ref.gff"

awk '$1 != "chrY" && $1 !~ /_alt$/ && $1 != "chrM"' $in_fn_2 > \
    "${out_dir}/ref.noY+M+alt.gff"
awk '$1 != "chrX" && $1 !~ /_alt$/ && $1 != "chrM"' $in_fn_2 > \
    "${out_dir}/ref.noX+M+alt.gff"