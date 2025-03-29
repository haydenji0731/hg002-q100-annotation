from tidy.utils import *

def load_tx_spans(fn) -> dict:
    hdr = ['ctg', 'src', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(fn, sep='\t', header=None, comment='#')
    df.columns = hdr
    tx_spans = dict()
    for _, row in df.iterrows():
        gfeat = gFeat(row, att_sep='=', schema=['ID', 'Parent', 'transcript_biotype'])
        pid = gfeat.parent
        tid = gfeat.fid
        if pid not in tx_spans:
            tx_spans[pid] = []
        tx_spans[pid].append((tid, gfeat.start, gfeat.end))
    return tx_spans
            
def qc_gene_spans(gan, tx_spans) -> None:
    for gid in gan.genes:
        gene = gan.genes[gid]
        min_tstart = None
        max_tend = None
        for tid in gene.children:
            tx = gan.txes[tid]
            if min_tstart is None:
                min_tstart = tx.start
            else:
                min_tstart = min(tx.start, min_tstart)
            if max_tend is None:
                max_tend = tx.end
            else:
                max_tend = max(tx.end, max_tend)
        gene.start = min_tstart
        gene.end = max_tend

        if gid in tx_spans:
            gstart = gene.start
            gend = gene.end
            for _, tstart, tend in tx_spans[gid]:
                if tstart < gstart:
                    gstart = tstart
                if tend > gend:
                    gend = tend
            gene.start = gstart
            gene.end = gend

def main(args):
    if not is_pickled(args.in_file):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: gan = pickle.load(f)

    tx_spans = load_tx_spans(args.input_2)
    qc_gene_spans(gan, tx_spans)

    s = gan.to_str(args.fmt)
    qced_fn = os.path.join(args.out_dir, f'qced.{args.fmt}')
    with open(qced_fn, 'w') as fh: fh.write(s)

    out_fn = os.path.join(args.out_dir, 'final.gff')
    cmd = f'cat {qced_fn} {args.input_2} > {out_fn}'
    print(cmd)
    call(cmd, shell=True)
    
    out_fn_2 = os.path.join(args.out_dir, 'final.sorted.gff')
    cmd = f'gffread -O -F --keep-exon-attrs {out_fn} > {out_fn_2}'
    print(cmd)
    call(cmd, shell=True)

    out_fn_3 = os.path.join(args.out_dir, 'final.sorted.fmted.gff')
    cmd = f"sed 's/geneID=/Parent=/g' {out_fn_2} > {out_fn_3}"
    print(cmd)
    call(cmd, shell=True)