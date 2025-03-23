from tidy.utils import *

def create_gfeat(args):
    row, att_sep, schema = args
    return gFeat(row, att_sep, schema)

def load_gan(args):
    if args.fmt.lower() == 'gff':
        att_sep = "="
    else:
        att_sep = " "
    hdr = ['ctg', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(args.in_file, sep='\t', header=None)
    df.columns = hdr

    cds_df = df[df['type'] == 'CDS']
    cds_feats = [gFeat(row, att_sep, args.schema) for _, row in cds_df.iterrows()]
    cds_parents = set([x.parent for x in cds_feats])
    print(tmessage(f"detected CDS records for {len(cds_parents)} transcripts", Mtype.PROG))

    exon_df = df[df['type'] == 'exon']
    exon_feats = [gFeat(row, att_sep, args.schema) for _, row in exon_df.iterrows()]
    pc_exon_feats = [x for x in exon_feats if x.parent in cds_parents] # only keep coding exons
    pc_exon_parents = set([x.parent for x in pc_exon_feats])
    assert len(cds_parents) == len(pc_exon_parents)
    print(tmessage(f"coding exon records loaded", Mtype.PROG))

    tx_df = df[df['type'] == 'transcript']
    tx_feats = [gFeat(row, att_sep, args.schema) for _, row in tx_df.iterrows()]
    pc_tx_feats = [x for x in tx_feats if x.fid in cds_parents] # only keep coding txes
    assert len(cds_parents) == len(pc_tx_feats)
    pc_tx_parents = set([x.parent for x in pc_tx_feats])
    print(tmessage(f"coding transcript records loaded", Mtype.PROG))

    gene_df = df[df['type'] == 'gene']
    gene_feats = [gFeat(row, att_sep, args.schema) for _, row in gene_df.iterrows()]
    pc_gene_feats = [x for x in gene_feats if x.fid in pc_tx_parents]
    assert len(pc_gene_feats) == len(pc_tx_parents)
    print(tmessage(f"coding gene records loaded", Mtype.PROG))

    if args.btype: # check biotypes if requested
        tx_btype_tbl = dict()
        for x in pc_tx_feats:
            btype = x.ftype
            if btype not in tx_btype_tbl: tx_btype_tbl[btype] = 1
            else: tx_btype_tbl[btype] += 1
        print("### biotype decomposition for coding transcripts ###")
        for k, v in tx_btype_tbl.items(): print(f'{k}\t{v}')

    return cds_feats, pc_exon_feats, pc_tx_feats, pc_gene_feats

def add_exon_ids(e_chain, parent):
    ctr = 1
    for e in e_chain:
        e.att_tbl['exon_number'] = ctr
        e.att_tbl['ID'] = f'{parent}-exon-{ctr}'
        ctr += 1

def add_cds_enumber(c_chain, e_chain):
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

def add_cds_ids(c_chain, parent):
    ctr = 1
    for c in c_chain:
        c.att_tbl['ID'] = f'{parent}-cds-{ctr}'
        ctr += 1

def build_gan(cdses, exons, txes, genes, add_id) -> gAn:
    print(tmessage(f"building gan", Mtype.PROG))
    gan = gAn(genes, txes)

    # build exon chains
    for x in exons:
        tid = x.parent
        tx = gan.txes[tid]
        if not tx.is_chain_init():
            tx.init_chain()
        tx.append2chain(0, x)
    
    # build cds chains
    for x in cdses:
        tid = x.parent
        tx = gan.txes[tid]
        assert tx.is_chain_init() # sanity check
        tx.append2chain(1, x)
    
    # add tx children
    for x in txes:
        gid = x.parent
        gan.genes[gid].add_child(x.fid) # just store tx ids
        # also sort chains
        gan.txes[x.fid].sort_chain(0)
        gan.txes[x.fid].sort_chain(1)
        if add_id:
            # add exon and cds IDs
            add_exon_ids(x.chains[0], x.fid)
            add_cds_ids(x.chains[1], x.fid)
        add_cds_enumber(x.chains[1], x.chains[0])
    
    return gan

def write_rest(in_fn, out_fn, gan, fmt):
    out_fh = open(out_fn, 'w')
    with open(in_fn, 'r') as in_fh:
        for ln in in_fh:
            if ln[0] == '#': continue
            gline = gLine(ln.strip(), fmt)
            if gline.feature == 'gene':
                gid = gline.attributes['ID']
                if gid in gan.genes: continue
                out_fh.write(ln)
            elif gline.feature == 'transcript':
                tid = gline.attributes['ID']
                if tid in gan.txes: continue
                out_fh.write(ln)
            else:
                pid = gline.attributes['Parent']
                if pid in gan.txes: continue
                out_fh.write(ln)
    out_fh.close()
            
def main(args):

    cdses, exons, txes, genes = load_gan(args)
    gan = build_gan(cdses, exons, txes, genes, args.add_id)

    s = gan.to_str(args.fmt)
    with open(os.path.join(args.out_dir, f'in.coding.{args.fmt}'), 'w') as fh: fh.write(s)

    save_pth = os.path.join(args.out_dir, 'in.coding.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)

    write_rest(args.in_file, os.path.join(args.out_dir, f'rest.{args.fmt}'), gan, args.fmt.lower())