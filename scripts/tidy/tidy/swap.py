from tidy.utils import *

def load_map(fn) -> dict:
    df = pd.read_csv(fn)
    id_mapping = dict()
    for _, row in df.iterrows():
        id_mapping[row['id_2']] = row['id_1']
    return id_mapping

def build_oid2tid_tbl(gan) -> dict:
    oid2tid_tbl = dict()
    for tid in gan.txes:
        tx = gan.txes[tid]
        oid = tx.att_tbl['origin_ID']
        oid2tid_tbl[oid] = tid
    return oid2tid_tbl
        
def swap_arecs(gan, fn, fmt, id_mapping, oid2tid_tbl, split_tid, tx_ftype, note) -> int:
    to_swap_info = dict()
    load = False
    with open(fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            gline = gLine(ln.strip(), fmt)
            if gline.feature == tx_ftype:
                load = False
                assert 'ID' in gline.attributes
                tid = gline.attributes['ID']
                if split_tid:
                    tid = tid.split('-')[1]
                if tid in id_mapping:
                    to_swap_info[tid] = [(gline.start, gline.end, gline.strand), [], []]
                    load = True
            elif gline.feature == 'exon':
                if load:
                    to_swap_info[tid][1].append((gline.start, gline.end, gline.strand))
            elif gline.feature == 'CDS':
                if load:
                    to_swap_info[tid][2].append((gline.start, gline.end, gline.strand))
            else:
                load = False
    for tid in id_mapping:
        if tid not in to_swap_info:
            print(tmessage(f"{tid} not found in second input", Mtype.WARN))
    for o_tid in to_swap_info:
        tstart, tend, strand = to_swap_info[o_tid][0]
        exons = to_swap_info[o_tid][1]
        cdses = to_swap_info[o_tid][2]

        oid = id_mapping[o_tid]
        if oid not in oid2tid_tbl:
            print(tmessage(f"{oid} not found in first input", Mtype.WARN))
            continue
        tid = oid2tid_tbl[oid]
        tx = gan.txes[tid]
        tx.start = tstart
        tx.end = tend
        assert tx.strand == strand # sanity check

        new_exon_chain = []
        exon_template = tx.chains[0][0]
        for e_start, e_end, e_strand in exons:
            exon = copy.deepcopy(exon_template)
            exon.start = e_start
            exon.end = e_end
            assert exon.strand == e_strand
            new_exon_chain.append(exon)
        tx.set_chain(0, new_exon_chain)
        tx.sort_chain(0)

        new_cds_chain = []
        cds_template = tx.chains[1][0]
        for c_start, c_end, c_strand in cdses:
            cds = copy.deepcopy(cds_template)
            cds.start = c_start
            cds.end = c_end
            assert cds.strand == c_strand
            new_cds_chain.append(cds)
        tx.set_chain(1, new_cds_chain)
        tx.sort_chain(1)
        tx.att_tbl['manual'] = 'True'
        if note: tx.att_tbl['note'] = note
    
    return len(to_swap_info)

def qc_gan(gan) -> None:
    for tid in gan.txes:
        tx = gan.txes[tid]
        e_chain = tx.chains[0]
        c_chain = tx.chains[1]

        ctr = 1
        for e in e_chain:
            e.att_tbl['exon_number'] = ctr
            if "ID" in e.att_tbl:
                e.att_tbl['ID'] = f'{tid}-exon-{ctr}'
            ctr += 1

        starting_cds = c_chain[0]
        starting_ei = -1
        for i, e in enumerate(e_chain):
            if starting_cds.start >= e.start and starting_cds.end <= e.end:
                starting_ei = i
                break
        assert starting_ei != -1

        starting_ei += 1
        ctr = 1
        for c in c_chain:
            assert starting_ei <= len(e_chain) + 1
            c.att_tbl['exon_number'] = starting_ei
            starting_ei += 1
            if "ID" in c.att_tbl:
                c.att_tbl['ID'] = f'{tid}-cds-{ctr}'
            ctr += 1

        tstart = tx.start
        tend = tx.end
        gene = gan.genes[tx.parent]
        if tstart < gene.start:
            gene.start = tstart
        if tend < gene.end:
            gene.end = tend

def main(args):
    if not is_pickled(args.in_file):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: gan = pickle.load(f)
    oid2tid_tbl = build_oid2tid_tbl(gan)

    id_mapping = load_map(args.map_file)
    if args.refseq:
        split_tid = True
        tx_ftype = 'mRNA'
    else:
        split_tid = False
        tx_ftype = 'transcript'
    n_swapped = swap_arecs(gan, args.input_2, args.fmt.lower(), id_mapping, \
                        oid2tid_tbl, split_tid, tx_ftype, args.note)
    print(tmessage(f"swapped records for {n_swapped} transcripts", Mtype.PROG))
    qc_gan(gan)

    print(tmessage(f"saving results", Mtype.PROG))
    s = gan.to_str(args.fmt)
    with open(os.path.join(args.out_dir, f'swapped.{args.fmt}'), 'w') as fh: fh.write(s)
    save_pth = os.path.join(args.out_dir, 'swapped.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)