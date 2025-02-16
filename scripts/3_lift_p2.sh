#!/usr/bin/env bash

set -x
# prep rDNA reference & run 2nd liftoff

chroms_fn=$1
tgt_genome=$2
out_dir=$3
prefix=$4

ref_fn="/ccb/salz4-3/hji20/hg002-q100-annotation/data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"
ref_genome="/ccb/salz4-3/hji20/hg002-q100-annotation/data/chm13v2.0.fa"
spec_fn="/ccb/salz4-3/hji20/hg002-q100-annotation/data/id_spec.2.csv"
features_fn="/ccb/salz4-3/hji20/hg002-q100-annotation/data/features.2.txt"
threads=24

mkdir -p $out_dir
./make_rDNA_ref.py $ref_fn gff "${out_dir}/rDNA_units.gff" gff # no need to repeat this again
rdna_fn="${out_dir}/rDNA_units.gff"
./make_gffutils_db.py $rdna_fn "${rdna_fn}_db" $spec_fn
liftoff -db "${rdna_fn}_db" -chroms $chroms_fn -p $threads -f $features_fn \
    -o "${out_dir}/${prefix}.gff" -u "${out_dir}/unmapped.txt" \
    -copies -sc 0.95 -exclude_partial -polish -dir "${out_dir}/tmp" \
    $tgt_genome $ref_genome