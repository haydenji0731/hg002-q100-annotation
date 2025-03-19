from utils import *

def create_gfeat(args):
    row, att_sep, schema = args
    return gFeat(row, att_sep, schema)

def load_gan(args):
    if args.ext.lower() == 'gff':
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

def build_gan(cdses, exons, txes, genes) -> gAn:
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
    
    return gan
    
def main(args):

    cdses, exons, txes, genes = load_gan(args)
    gan = build_gan(cdses, exons, txes, genes)
    save_pth = os.path.join(args.out_dir, 'in.coding.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)