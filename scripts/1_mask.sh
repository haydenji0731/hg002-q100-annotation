#!/usr/bin/env bash

set -x 

# mask rDNA regions of the target genome

ctg_fn=$1
vdj_fn=$2
out_dir=$3
prefix=$4

full_genome="/ccb/salz4-3/hji20/hg002-q100-annotation/data/hg002v1.1.fasta"
rdna_fn="/ccb/salz4-3/hji20/hg002-q100-annotation/data/rDNA.fa"
similarity=96
length=1000
threads=36

samtools faidx -r $ctg_fn $full_genome > "${out_dir}/${prefix}.fa"
samtools faidx "${out_dir}/${prefix}.fa"

nucmer --maxmatch -t $threads -l 31 -c 100 "${out_dir}/${prefix}.fa" $rdna_fn -p "${out_dir}/${prefix}"
delta-filter -1 "${out_dir}/${prefix}.delta" > "${out_dir}/${prefix}.1to1.delta"
show-coords -lcHr "${out_dir}/${prefix}.1to1.delta" > "${out_dir}/${prefix}.1to1.coords"
show-coords -lcHr "${out_dir}/${prefix}.delta" > "${out_dir}/${prefix}.coords"

awk -v s="${similarity}" -v l="${length}" '{if($10 > s && $7 > l) {print}}' \
   < "${out_dir}/${prefix}.coords" > "${out_dir}/${prefix}.filtered.coords"
awk '{print $18 "\t" $1-1 "\t" $2}' < "${out_dir}/${prefix}.filtered.coords" > "${out_dir}/${prefix}.filtered.bed"

cat $vdj_fn "${out_dir}/${prefix}.filtered.bed" > "${out_dir}/${prefix}.mask_regions.bed"

bedtools maskfasta -fi "${out_dir}/${prefix}.fa" -bed "${out_dir}/${prefix}.mask_regions.bed" -fo "${out_dir}/${prefix}.masked.fasta"
samtools faidx "${out_dir}/${prefix}.masked.fasta"