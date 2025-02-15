#!/usr/bin/env bash

set -x

ref_fn=$1
chroms_fn=$2
tgt_genome=$3
out_dir=$4
prefix=$5

ref_genome="/ccb/salz4-3/hji20/hg002-q100-annotation/data/chm13v2.0.fa"
spec_fn="/ccb/salz4-3/hji20/hg002-q100-annotation/data/id_spec.csv"
threads=24

mkdir -p $out_dir
./make_gffutils_db.py $ref_fn "${ref_fn}_db" $spec_fn

# mm2 -N 100 set as default
liftoff -db "${ref_fn}_db" -chroms $chroms_fn \
    -p $threads -o "${out_dir}/${prefix}.gff" -u "${out_dir}/unmapped.txt" \
    -copies -sc 0.95 -exclude_partial -polish -dir "${out_dir}/tmp" \
    $tgt_genome $ref_genome