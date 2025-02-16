#!/usr/bin/env bash

# remove rDNA unit fake roots and merge p1 and p2 outputs

set -x

p1_fn=$1
p2_fn=$2
out_dir=$3
prefix=$4

mkdir -p $out_dir
./remove_fake_root.py $p2_fn gff "${out_dir}/${prefix}.p2.fmted.gff" gff
p2_clean_fn="${out_dir}/${prefix}.p2.fmted.gff"

# check overlapping features
awk -v OFS='\t' '{ print $1, $4, $5 }' $p1_fn | sed '/^#/d' > "${out_dir}/${prefix}.p1.bed"
awk -v OFS='\t' '{ print $1, $4, $5 }' $p2_clean_fn | sed '/^#/d' > "${out_dir}/${prefix}.p2.bed"
bedtools intersect -wa -wb -a "${out_dir}/${prefix}.p1.bed" \
    -b "${out_dir}/${prefix}.p2.bed" > "${out_dir}/${prefix}.p1_2.intersect.bed"
awk -v OFS='\t' '{ print $4, $5, $6 }' "${out_dir}/${prefix}.p1_2.intersect.bed" | sort | uniq > "${out_dir}/${prefix}.p2.regions.bed"

# actual merge between gff files
cat $p1_fn $p2_clean_fn > "${out_dir}/${prefix}.final.gff"
gffread -O -F --keep-exon-attrs "${out_dir}/${prefix}.final.gff" > "${out_dir}/${prefix}.final.sorted.gff"