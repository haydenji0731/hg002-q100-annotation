#!/usr/bin/env python

# remove vdj segments, rDNA, and additional gene copies

from commons import *

# filters on biotypes
def sift(in_fn, format, gene_filters, tx_filters):
    remove_genes = []
    remove_txes = []
    with open(in_fn, 'r') as in_fh:
        for ln in in_fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, format)

            if ln_o.feature == 'gene':
                fid = ln_o.attributes['ID']
                # assert 'gene_biotype' in ln_o.attributes # sanity check; success
                if 'extra_copy_number' in ln_o.attributes and \
                    ln_o.attributes['extra_copy_number'] != '0':
                    remove_genes.append(fid)
                elif 'gene_biotype' in ln_o.attributes and \
                    ln_o.attributes['gene_biotype'] in gene_filters:
                    if ln_o.attributes['gene_biotype'] == 'rRNA' and 'RNA5S' in ln_o.attributes['ID']:
                        continue

                    remove_genes.append(fid)
            elif ln_o.feature == 'transcript':
                fid = ln_o.attributes['ID']
                if 'extra_copy_number' in ln_o.attributes and \
                    ln_o.attributes['extra_copy_number'] != '0':
                    remove_txes.append(fid)
                elif 'transcript_biotype' in ln_o.attributes and \
                    ln_o.attributes['transcript_biotype'] in tx_filters:
                    if ln_o.attributes['transcript_biotype'] == 'rRNA' and 'RNA5S' in ln_o.attributes['ID']:
                        continue
                    remove_txes.append(fid)
    return remove_genes, remove_txes

def write_filtered(in_fn, in_format, out_fn, out_format, remove_genes, remove_txes) -> dict:
    gene_fields = ['ID', 'gene_name', 'gene_biotype']
    tx_fields = ['ID', 'Parent', 'transcript_biotype']
    misc_fields = ['Parent', 'exon_number']
    out_fh = open(out_fn, 'w')
    ctr = 1
    prev_ctg = None
    old2new_gene = dict()
    isoform_ctr = dict()
    old2new_tx = dict()
    tsle = dict()
    with open(in_fn, 'r') as in_fh:
        for ln in in_fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, in_format)
            if ln_o.feature == 'gene':
                if ln_o.attributes['ID'] in remove_genes: continue
                if ln_o.ctg != prev_ctg:
                    prev_ctg = ln_o.ctg
                    ctr = 1
                else:
                    ctr += 1
                temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in gene_fields}
                temp['origin_ID'] = temp['ID']
                temp['ID'] = f'REF.{ln_o.ctg}.{ctr}'
                old2new_gene[temp['origin_ID']] = temp['ID']
                ln_o.attributes = temp
                out_fh.write(ln_o.to_gStr(out_format))
            elif ln_o.feature == 'transcript':
                if ln_o.attributes['ID'] in remove_txes: continue
                assert ln_o.attributes['Parent'] not in remove_genes
                temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in tx_fields}
                temp['origin_ID'] = temp['ID']
                pid = old2new_gene[temp['Parent']]
                if pid not in isoform_ctr:
                    isoform_ctr[pid] = 1
                else:
                    isoform_ctr[pid] += 1
                temp['ID'] = f'{pid}.{isoform_ctr[pid]}'
                old2new_tx[temp['origin_ID']] = temp['ID']
                temp['Parent'] = pid
                ln_o.attributes = temp
                out_fh.write(ln_o.to_gStr(out_format))
            elif ln_o.feature == 'exon':
                if ln_o.attributes['Parent'] in remove_txes: continue
                temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in misc_fields}
                pid = old2new_tx[temp['Parent']]
                temp['Parent'] = pid
                ln_o.attributes = temp
                out_fh.write(ln_o.to_gStr(out_format))
            elif ln_o.feature == 'CDS':
                if ln_o.attributes['Parent'] in remove_txes: continue
                temp = {k: ln_o.attributes[k] if k in ln_o.attributes else 'NA' for k in misc_fields}
                if 'transl_except' in ln_o.attributes:
                    tsle_temp = ln_o.attributes['transl_except'].split('%')
                    assert len(tsle_temp) == 2
                    if 'complement' in tsle_temp[0]:
                        pos = [x.replace(')', '') for x in tsle_temp[0].split('(')[2].split('..')]
                    else:
                        pos = [x for x in tsle_temp[0].split(':')[1].split('..')]
                    excpt = tsle_temp[1].replace('2C', '').replace(')', '')
                    assert excpt[:2] == 'aa' # sanity check
                    if temp['Parent'] not in tsle:
                        tsle[temp['Parent']] = (pos, excpt)
                pid = old2new_tx[temp['Parent']]
                temp['Parent'] = pid
                ln_o.attributes = temp
                out_fh.write(ln_o.to_gStr(out_format))
    out_fh.close()
    return tsle

def main(in_fn, out_dir, in_format, out_format):
    remove_genes = []
    remove_txes = []
    gene_filters = ['V_segment', 'C_segment', 'J_segment', 'rRNA']
    tx_filters = ['V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'rRNA']
    remove_genes, remove_txes = sift(in_fn, in_format, gene_filters, tx_filters)
    print(tmessage(f'{len(remove_genes)} genes to remove', Mtype.PROG))
    print(tmessage(f'{len(remove_txes)} transcripts to remove', Mtype.PROG))
    tsle = write_filtered(in_fn, in_format, os.path.join(out_dir, 'ref.gff'), out_format, remove_genes, remove_txes)
    with open(os.path.join(out_dir, 'ref.tsl.except'), 'w') as fh:
        fh.write('transcript_id,pos,exception\n')
        for tid in tsle:
            fh.write(f'{tid},{":".join(tsle[tid][0])},{tsle[tid][1]}\n')

    print(tmessage(f'finished cleaning', Mtype.PROG))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])