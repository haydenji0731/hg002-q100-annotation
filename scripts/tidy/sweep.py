from utils import *

def extend_upstream(s, cds_start, w) -> int:
    wd_end = max(0, cds_start - w) # respect 5' end
    for i in range(cds_start, wd_end, -3):
        if s[i-3:i] in START_CODONS:
            return i-3
        if s[i-3:i] in STOP_CODONS:
            return cds_start
    return cds_start

def extend_downstream(s, cds_end, w) -> int:
    wd_end = min(len(s), cds_end + w) # respect 3' end
    for i in range(cds_end, wd_end, 3):
        if s[i:i+3] in STOP_CODONS:
            return i+3
    return cds_end

def check_cds_integrity(cds_chain):
    is_valid = True
    for cds in cds_chain:
        if cds.end < cds.start:
            is_valid = False
            break
    return is_valid

# adopted from shorten_chain fn in trim.py
def trim_chain(in_chain, strand, l, side):
    out_chain = copy.deepcopy(in_chain)
    remaining = l
    if strand == '+':
        if side == 0: # start
            while remaining > 0:
                curr_clen = out_chain[0].end - out_chain[0].start + 1
                if remaining >= curr_clen:
                    out_chain = out_chain[1:]
                    remaining -= curr_clen
                else:
                    out_chain[0].start += remaining
                    remaining -= remaining # sets to zero
        else:
            while remaining > 0:
                curr_clen = out_chain[-1].end - out_chain[-1].start + 1
                if remaining >= curr_clen:
                    out_chain = out_chain[:-1]
                    remaining -= curr_clen
                else:
                    out_chain[-1].end -= remaining
                    remaining -= remaining
    else:
        if side == 0: # start
            while remaining > 0:
                curr_clen = out_chain[0].end - out_chain[0].start + 1
                if remaining >= curr_clen:
                    out_chain = out_chain[1:]
                    remaining -= curr_clen
                else:
                    out_chain[0].end -= remaining
                    remaining -= remaining
        else:
            while remaining > 0:
                curr_clen = out_chain[-1].end - out_chain[-1].start + 1
                if remaining >= curr_clen:
                    out_chain = out_chain[:-1]
                    remaining -= curr_clen
                else:
                    out_chain[-1].start += remaining
                    remaining -= remaining
    assert check_cds_integrity(out_chain)
    return out_chain

def define_cds(start, end, tx, n) -> None:
    assert tx.sorted[0] # check if exon chain is sorted
    temp_chain = trim_chain(copy.deepcopy(tx.chains[0]), tx.strand, start, 0)
    new_cds_chain = trim_chain(temp_chain, tx.strand, n - end, 1)
    tx.set_chain(1, new_cds_chain)

def find_alt_start_stop(tx, cds_coords, nseq, w):
    if tx.att_tbl['valid_ORF'] == 'True': return False, None
    else:
        if 'inframe_stop_codon' in tx.att_tbl: return False, None
        if 'missing_start_codon' in tx.att_tbl:
            start_pos = extend_upstream(nseq, cds_coords[0], w)
        else:
            start_pos = cds_coords[0]
        if 'missing_stop_codon' in tx.att_tbl:
            end_pos = extend_downstream(nseq, cds_coords[1], w)
        else:
            end_pos = cds_coords[1]
        if start_pos == cds_coords[0] and end_pos == cds_coords[1]: return False, None
        define_cds(start_pos, end_pos, tx, len(nseq))
        if nseq[start_pos:start_pos + 3] in START_CODONS and nseq[end_pos - 3:end_pos] in STOP_CODONS:
            tx.att_tbl['valid_ORF'] = 'True'
        tx.att_tbl['manual'] = 'True'
    return True, (tx.fid, start_pos, end_pos, cds_coords[0], cds_coords[1])

# TODO: implement
def is_orf_contained(nseq, cds_coords, ref_pseq):
    cds_start, cds_end = cds_coords
    assert (cds_end - cds_start) % 3 == 0
    pseq = Seq.translate(nseq[cds_start:cds_end])
    ref_pseq_wstop = ref_pseq + '*'
    if ref_pseq_wstop == pseq:
        return False, None
    if ref_pseq in pseq:
        start = pseq.find(ref_pseq)
        end = start + len(ref_pseq)
        return True, (start, end)
    return False, None

# TODO: implement
def find_alt_orf():
    return

def load_cds_coords(fn) -> dict:
    cds_coords_tbl = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            clean_ln = ln.strip()
            if clean_ln[0] == '>':
                temp = clean_ln.split(" ")
                if len(temp) == 2:
                    cds_info = temp[1]
                    coords = [int(x) for x in cds_info.split('=')[1].strip().split('-')]
                    tid = temp[0].replace('>', '')
                    coords[0] -= 1 # 1-based to 0-based, half-open system
                    cds_coords_tbl[tid] = coords
    return cds_coords_tbl

def main(args):
    if not is_pickled(args.in_file):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    if args.window % 3 != 0:
        print(tmessage("extend window size must be divisible by 3", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: gan = pickle.load(f)

    if args.slip_except:
        print(tmessage("loading ribosomal slippage exceptions", Mtype.PROG))
        if not os.path.exists(args.slip_except):
            print(tmessage("slippage exception file not found", Mtype.ERR))
            sys.exit(-1)
        rslippage_tids = load_lst(args.slip_except)
    else:
        rslippage_tids = None

    nfa = pyfastx.Fasta(args.nucl_file)
    cds_coords_tbl = load_cds_coords(args.nucl_file)
    rfa = pyfastx.Fasta(args.ref_file)
    # assert len(cds_coords_tbl) == len(gan.txes) # sanity check

    print(tmessage("finding alternative start and/or stop codons", Mtype.PROG))
    new_cdses = []
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        is_new_cds, coords = find_alt_start_stop(tx, cds_coords_tbl[tid], nfa[tid].seq, args.window)
        if is_new_cds:
            new_cdses.append(coords) # tid, new_start, new_end, old_start, old_end
            cds_coords_tbl[tid] = (coords[1], coords[2])
    write_tup_lst(new_cdses, os.path.join(args.out_dir, 'start_stop_adj.csv'))
    print(tmessage(f"alternative start and/or stop found for {len(new_cdses)} transcripts", Mtype.PROG))

    print(tmessage("checking reference orf containments", Mtype.PROG))
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        if tx.att_tbl['origin_ID'] in rslippage_tids: continue
        contained, coords = is_orf_contained(nfa[tid].seq, cds_coords_tbl[tid], rfa[tx.att_tbl['origin_ID']].seq)
        if contained:
            # TODO: handle containment by defining a new CDS like above
            continue
    
    

        

