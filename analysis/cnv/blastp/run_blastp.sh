#!/usr/bin/env bash

set -x

blastp="/ccb/salz4-3/hji20/hg002-q100-annotation/tools/ncbi-blast-2.16.0+/bin/blastp"
# NOTE: MUST review this before running
out_dir="/ccb/salz4-3/hji20/hg002-q100-annotation/results/analysis/blastp/pat_loff"

ref_db=$1
shift

for qry_fn in "$@"; do
    echo "processing $qry_fn"
    qry=$(basename "$qry_fn" .fa)
    $blastp -query $qry_fn \
            -db $ref_db \
            -out "${out_dir}/${qry}.blastp.out" \
            -outfmt "6 std slen qlen"
done