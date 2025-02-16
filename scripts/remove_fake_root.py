#!/usr/bin/env python

# removes fake roots from rDNA liftoff results

from commons import *

def main(in_fn, in_format, out_fn, out_format):
    out_fh = open(out_fn, 'w')
    ctr = 1
    old2new_tx = dict()
    with open(in_fn, 'r') as in_fh:
        for ln in in_fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, in_format)
            if ln_o.feature == 'unit':
                continue
            elif ln_o.feature == 'transcript':
                gid = f'RDNA.{ctr}'
                gene_o = copy_o(ln_o, gid, in_format, 'gene')
                gene_o.attributes['gene_name'] = ln_o.attributes['gene_name']
                gene_o.attributes['origin_ID'] = ln_o.attributes['origin_ID'].replace('T', 'G')
                gene_o.attributes['unit_ID'] = ln_o.attributes['Parent']
                gene_o.attributes['gene_biotype'] = "rRNA"
                gene_o.attributes['extra_copy_number'] = ln_o.attributes['extra_copy_number']
                out_fh.write(gene_o.to_gStr(out_format))

                old2new_tx[ln_o.attributes['ID']] = f'{gid}.1'
                ln_o.attributes['ID'] = f'{gid}.1'
                ln_o.attributes['Parent'] = gid
                out_fh.write(ln_o.to_gStr(out_format))
                ctr += 1
            else: # exons
                ln_o.attributes.pop('ID') # make it consistent with liftoff's output
                ln_o.attributes['Parent'] = old2new_tx[ln_o.attributes['Parent']]
                out_fh.write(ln_o.to_gStr(out_format))
    out_fh.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])