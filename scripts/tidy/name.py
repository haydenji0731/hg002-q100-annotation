from utils import *

def main(args):
    prev_ctg = None
    isoform_ctr = dict()
    old2new_gene = dict()
    old2new_tx = dict()
    out_fh = open(args.out_file, 'w')
    file_format = args.ext.lower()
    with open(args.in_file, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, file_format)
            if ln_o.feature == "gene":
                if ln_o.ctg != prev_ctg:
                    prev_ctg = ln_o.ctg
                    gene_ctr = 1
                else:
                    gene_ctr += 1
                old_gid = ln_o.attributes['ID']
                new_gid = f'{args.prefix}_{ln_o.ctg[:args.char_limit]}_{gene_ctr}'
                ln_o.attributes['ID'] = new_gid
                ln_o.attributes['copy_num_ID'] = f'{new_gid}_{ln_o.attributes['extra_copy_number']}'
                assert old_gid not in isoform_ctr
                isoform_ctr[old_gid] = 1
                old2new_gene[old_gid] = new_gid
                out_fh.write(ln_o.to_gStr(file_format))
            elif ln_o.feature == "transcript":
                pid = ln_o.attributes['Parent']
                assert pid in isoform_ctr
                old_tid = ln_o.attributes['ID']
                new_tid = f'{old2new_gene[pid]}.{isoform_ctr[pid]}'
                ln_o.attributes['ID'] = new_tid
                ln_o.attributes['Parent'] = old2new_gene[pid]
                old2new_tx[old_tid] = new_tid
                isoform_ctr[pid] += 1
                out_fh.write(ln_o.to_gStr(file_format))
            elif ln_o.feature == 'exon' or ln_o.feature == 'CDS':
                pid = ln_o.attributes['Parent']
                assert pid in old2new_tx
                ln_o.attributes['Parent'] = f'{old2new_tx[pid]}'
                out_fh.write(ln_o.to_gStr(file_format))
            else:
                print(ln_o.feature)
                print(tmessage(f'invalid feature', Mtype.ERR))
                sys.exit(-1)
    out_fh.close()