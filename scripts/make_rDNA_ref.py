#!/usr/bin/env python

# isolate rDNA reference from CHM13 annotation

from commons import *

def isolate_rdna_units(fn, format) -> dict:
    rdna_subtypes = ["RNA18S", "RNA5-8S", "RNA28S"]
    gene_fields = ['ID', 'gene_name', 'gene_biotype']
    unit_ctr = 0
    unit_d = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, format)
            if ln_o.feature == 'gene':
                if ln_o.attributes['gene_biotype'] == 'rRNA':
                    if 'RNA45S' in ln_o.attributes['gene_name']:
                        unit_ctr += 1
                        el_ctr = 1
                        unit_id = f'RDNA.{unit_ctr}' # ctg info omitted
                        unit_o = copy_o(ln_o, unit_id, format, "unit")
                        unit_d[unit_id] = (unit_o, [])
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
           for x in genes:
               x.feature = "transcript" # change feature type to transcript from gene
               x.attributes['origin_ID'] = x.attributes['origin_ID'].replace('G', 'T')
               x.attributes['transcript_biotype'] = x.attributes.pop('gene_biotype')
               x.attributes.pop("gene_name")
               fh.write(x.to_gStr(format))
               exon_o = copy_o(x, f'{x.attributes["ID"]}-exon-1', format, "exon")
               exon_o.attributes['Parent'] = x.attributes["ID"]
               exon_o.attributes['exon_number'] = 1
               fh.write(exon_o.to_gStr(format))

def main(in_fn, in_format, out_fn, out_format):
    unit_d = isolate_rdna_units(in_fn, in_format)
    print(tmessage(f'identified {len(unit_d)} full rDNA units', Mtype.PROG))
    write_rdna(unit_d, out_fn, out_format)
    print(tmessage(f'finished', Mtype.PROG))
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])