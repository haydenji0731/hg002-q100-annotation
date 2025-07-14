from tidy.utils import *

def add_cds_enumber(c_chain, e_chain) -> gAn:
    # find starting exon
    starting_cds = c_chain[0]
    starting_ei = -1
    for i, e in enumerate(e_chain):
        if starting_cds.start >= e.start and starting_cds.end <= e.end:
            starting_ei = i
            break
    assert starting_ei != -1
    starting_ei += 1 # switch to 1-based index
    for c in c_chain:
        assert starting_ei <= len(e_chain) + 1
        c.att_tbl['exon_number'] = starting_ei
        starting_ei += 1

def add_exon_number(e_chain) -> None:
    ctr = 1
    for x in e_chain:
        x.att_tbl['exon_number'] = ctr
        ctr += 1

def check_orf_validity(cseq):
    missing_start = True
    missing_stop = True
    inframe_stop = False
    l_violation = True
    start_codon = cseq[:3].upper()
    stop_codon = cseq[-3:].upper()
    if start_codon in START_CODONS:
        missing_start = False
    if stop_codon in STOP_CODONS:
        missing_stop = False
    if '*' in Seq.translate(cseq)[:-1]:
        inframe_stop = True
    if len(cseq) % 3 == 0:
        l_violation = False
    # print(f'{start_codon}\t{stop_codon}\t{missing_start}\t{missing_stop}')
    return missing_start, missing_stop, inframe_stop, l_violation

# load the gan to add to the existing
def load_and_build_gan(fn, fmt, schema, cfn):
    cfa = pyfastx.Fasta(cfn)
    orf_infos = dict()

    if fmt.lower() == 'gff':
        att_sep = '='
    else:
        att_sep = ' '
    
    hdr = ['ctg', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(fn, sep='\t', header=None, comment='#')
    df.columns = hdr

    genes = []
    txes = []
    exons = []
    cdses = []
    for _, row in df.iterrows():
        feat = gFeat(row, att_sep, schema)
        if row['type'] == 'gene':
            genes.append(feat)
        elif row['type'] == 'transcript':
            txes.append(feat)
        elif row['type'] == 'exon':
            exons.append(feat)
        elif row['type'] == 'CDS':
            cdses.append(feat)

    gan = gAn(genes, txes)

    for x in exons:
        tid = x.parent
        tx = gan.txes[tid]
        if not tx.is_chain_init(): tx.init_chain()
        tx.append2chain(0, x)
    
    for x in cdses:
        tid = x.parent
        tx = gan.txes[tid]
        if not tx.is_chain_init(): tx.init_chain()
        tx.append2chain(1, x)
    
    for x in txes:
        gid = x.parent
        gan.genes[gid].add_child(x.fid) # just store tx ids
        
        assert x.fid in cfa
        missing_start, missing_stop, inframe_stop, l_violation = check_orf_validity(cfa[x.fid].seq)
        orf_infos[x.fid] = (missing_start, missing_stop, inframe_stop, l_violation)

        # also sort chains
        gan.txes[x.fid].sort_chain(0)
        gan.txes[x.fid].sort_chain(1)
        add_exon_number(gan.txes[x.fid].chains[0])
        add_cds_enumber(gan.txes[x.fid].chains[1], gan.txes[x.fid].chains[0])
        
    return gan, orf_infos

def fmt(in_fn, out_fn, sub_gan, sub_oinfos, fmt, id_prefix, plc_holder) -> None:
    prev_ctg = None
    isoform_ctr = dict()
    old2new_gene = dict()
    old2new_tx = dict()

    out_fh = open(out_fn, 'w')

    with open(in_fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_o = gLine(ln, fmt)
            if ln_o.feature == 'gene':
                if ln_o.ctg != prev_ctg:
                    prev_ctg = ln_o.ctg
                    gene_ctr = 1
                else:
                    gene_ctr += 1
                old_gid = ln_o.attributes['ID']
                new_gid = f'{id_prefix}_{ln_o.ctg.lower().replace("x", "X").replace("y", "Y")}_{gene_ctr}'
                if old_gid in sub_gan.genes and ln_o.src == 'miniprot':
                    tid = list(sub_gan.genes[old_gid].children)[0]

                    temp = {
                        'ID': new_gid,
                        'gene_name': ln_o.attributes['gene_name'] if plc_holder in old_gid else ln_o.attributes['gene'],
                        'gene_biotype': 'protein_coding',
                        # next three attributes might need manual curation
                        'origin_ID': old_gid,
                        'valid_ORFs': 1 if True not in sub_oinfos[tid] else 0,
                        'extra_copy_number': 0,
                        'ref': ln_o.attributes['source'],
                    }

                    ln_o.attributes = temp
                else:
                    ln_o.attributes['ID'] = new_gid
                isoform_ctr[old_gid] = 1
                old2new_gene[old_gid] = new_gid
                out_fh.write(ln_o.to_gStr(fmt))
            elif ln_o.feature == 'transcript':
                pid = ln_o.attributes['Parent']
                assert pid in isoform_ctr
                old_tid = ln_o.attributes['ID']
                new_tid = f'{old2new_gene[pid]}.{isoform_ctr[pid]}'

                if old_tid in sub_gan.txes and ln_o.src == 'miniprot':
                    is_valid_orf = True if True not in sub_oinfos[old_tid] else False
                    missing_start, missing_stop, inframe_stop, l_violation = sub_oinfos[old_tid]
                    temp = {
                        'ID': new_tid,
                        'Parent': old2new_gene[pid],
                        'transcript_biotype': 'mRNA',
                        'origin_ID': old_tid,
                        'matches_ref_protein': 'NA',
                        'valid_ORF': is_valid_orf,
                        'extra_copy_number': 0,
                        'manual': 'True',
                        'no_utr': 'True',
                    }
                    if missing_start: temp['missing_start'] = 'True'
                    if missing_stop: temp['missing_stop'] = 'True'
                    if inframe_stop: temp['inframe_stop_codon'] = 'True'
                    if l_violation: temp['cds_not_div3'] = 'True'
                    
                    ln_o.attributes = temp
                else:
                    ln_o.attributes['ID'] = new_tid
                    
                old2new_tx[old_tid] = new_tid
                isoform_ctr[pid] += 1
                out_fh.write(ln_o.to_gStr(fmt))
            elif ln_o.feature == 'exon':
                pid = ln_o.attributes['Parent']
                assert pid in old2new_tx
                ln_o.attributes['Parent'] = f'{old2new_tx[pid]}'
                out_fh.write(ln_o.to_gStr(fmt))
            elif ln_o.feature == 'CDS':
                pid = ln_o.attributes['Parent']
                assert pid in old2new_tx
                ln_o.attributes['Parent'] = f'{old2new_tx[pid]}'
                if pid in sub_gan.txes and ln_o.src == 'miniprot':
                    ln_o.attributes['miniprot_ID'] = ln_o.attributes['miniprot_Identity']
                    del ln_o.attributes['miniprot_Identity']
                    if 'StopCodon' in ln_o.attributes: del ln_o.attributes['StopCodon']
                out_fh.write(ln_o.to_gStr(fmt))
    out_fh.close()
            
def main(args) -> None:
    print(tmessage(f'loading sub gan', Mtype.PROG))
    sub_gan, sub_orf_infos = load_and_build_gan(args.sub_file, args.fmt, args.schema, args.cds_file)
    print(tmessage(f'reformatting the input gan', Mtype.PROG))
    fmt(args.in_file, os.path.join(args.out_dir, 'fmted.gff'), sub_gan, sub_orf_infos, \
        args.fmt.lower(), args.id_prefix, args.plc_holder.lower())