#!/usr/bin/env python

# isolate rDNA reference from CHM13 annotation

from commons import *

def build_rdna_unit(gene_o, unit_id, format):
    unit_o = gLine(None, format)
    unit_o.ctg = gene_o.ctg
    unit_o.src = gene_o.src
    unit_o.start = gene_o.start
    unit_o.end = gene_o.end
    unit_o.strand = gene_o.strand
    unit_o.frame = gene_o.frame
    unit_o.score = gene_o.score
    unit_o.feature = "unit"
    unit_o.attributes["ID"] = unit_id
    return unit_o

def isolate_rdna_units(fn, format) -> dict:
    rdna_subtypes = ["RNA18S", "RNA5-8S", "RNA28S"]
    gene_fields = ['ID', 'gene_name', 'gene_biotype']
    unit_ctr = 1
    unit_d = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, format)
            if ln_o.feature == 'gene':
                if ln_o.attributes['gene_biotype'] == 'rRNA':
                    if 'RNA45S' in ln_o.attributes['gene_name']:
                        el_ctr = 1
                        unit_id = f'RDNA.{unit_ctr}' # ctg info omitted
                        unit_o = build_rdna_unit(ln_o, unit_id, format)
                        unit_d[unit_id] = (unit_o, [])
                        unit_ctr += 1
                        temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in gene_fields}
                        temp['Parent'] = unit_id
                        temp['origin_ID'] = temp['ID']
                        temp['ID'] = f'RDNA.{unit_ctr}.{el_ctr}'
                        ln_o.attributes = temp
                        unit_d[unit_id][1].append(ln_o)
                        el_ctr += 1
                    elif True in [x in ln_o.attributes['gene_name'] for x in rdna_subtypes]:
                        temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in gene_fields}
                        temp['Parent'] = unit_id
                        temp['origin_ID'] = temp['ID']
                        temp['ID'] = f'RDNA.{unit_ctr}.{el_ctr}'
                        ln_o.attributes = temp
                        unit_d[unit_id][1].append(ln_o)
                        el_ctr += 1
    # NOTE: gene, transcript, exon share start and end coordinates
    return unit_d

def write_rdna(unit_d, fn, format):    
    with open(fn, 'w') as fh:
       for unit_id in  unit_d:
           unit_o, genes = unit_d[unit_id]
           fh.write(unit_o.to_gStr(format))
           for gene_o in genes:
               fh.write(gene_o.to_gStr(format))

def main(in_fn, in_format, out_fn, out_format):
    unit_d = isolate_rdna_units(in_fn, in_format)
    write_rdna(unit_d, out_fn, out_format)
    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])