from tidy.utils import *
import orfipy_core

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

        if nseq[start_pos:start_pos + 3] in START_CODONS:
            tx.att_tbl.pop('missing_start_codon', None)
            valid_start = True
        else:
            valid_start = False
        if nseq[end_pos - 3:end_pos] in STOP_CODONS:
            tx.att_tbl.pop('missing_stop_codon', None)
            valid_stop = True
        else:
            valid_stop = False
        if valid_start and valid_stop:
            tx.att_tbl['valid_ORF'] = 'True'

        tx.att_tbl['manual'] = 'True'
    return True, (tx.fid, start_pos, end_pos, cds_coords[0], cds_coords[1])

def check_containment(tx, nseq, cds_coords, rseq):
    contained, coords = is_orf_contained(nseq, cds_coords, rseq)
    if contained:
        _, _, nstart, nend = coords
        define_cds(nstart, nend, tx, len(nseq))
        if nseq[nstart:nstart + 3] in START_CODONS:
            tx.att_tbl.pop('missing_start_codon', None)
            valid_start = True
        else:
            valid_start = False
        if nseq[nend - 3:nend] in STOP_CODONS:
            tx.att_tbl.pop('missing_stop_codon', None)
            valid_stop = True
        else:
            valid_stop = False
        if valid_start and valid_stop:
            tx.att_tbl['valid_ORF'] = 'True'
        else:
            tx.att_tbl['valid_ORF'] = 'False'
        tx.att_tbl['manual'] = 'True'
        return contained, (tx.fid, nstart, nend, cds_coords[0], cds_coords[1])
    return contained, None

def is_orf_contained(nseq, cds_coords, rseq):
    cds_start, cds_end = cds_coords
    assert (cds_end - cds_start) % 3 == 0
    pseq = Seq.translate(nseq[cds_start:cds_end])
    rseq_wstop = rseq + '*'
    if rseq_wstop == pseq or rseq == pseq:
        return False, None
    # if rseq in pseq, then it's simply a longer protein
    if rseq_wstop in pseq:
        pstart = pseq.find(rseq_wstop)
        pend = pstart + len(rseq_wstop)
        nstart = 3 * pstart + cds_start
        nend = 3 * pend + cds_start
        return True, (pstart, pend, nstart, nend)
    return False, None

# NOTE: useful when hg38 reference is provided
def has_alt_orf(nseq, rseq, cds_coords, ratio, min_match):
    rlen = len(rseq) * 3
    rseq_wstop = rseq + '*' # NOTE: it's not guaranteed that rseq ends with a valid stop
    for start, stop, _, _ in orfipy_core.orfs(nseq, minlen=ratio * rlen, \
                                                    maxlen=rlen, strand='f', \
                                                    starts=START_CODONS, stops=STOP_CODONS):
        if cds_coords[0] == start and cds_coords[1] == stop + 3: continue
        alt_cseq = nseq[start:stop + 3]
        tsl_alt_cseq = Seq.translate(alt_cseq)

        if tsl_alt_cseq == rseq_wstop or \
            (tsl_alt_cseq[:min_match] == rseq_wstop[:min_match] and \
            tsl_alt_cseq[-min_match:] == rseq_wstop[-min_match:]):
            return True, (start, stop + 3)
    return False, None

def find_alt_orf(tx, nseq, cds_coords, rseq, ratio=0.9, min_match=10):
    has_ao, coords = has_alt_orf(nseq, rseq, cds_coords, ratio, min_match)
    if has_ao:
        define_cds(coords[0], coords[1], tx, len(nseq))
        if nseq[coords[0]:coords[0] + 3] in START_CODONS:
            tx.att_tbl.pop('missing_start_codon', None)
            valid_start = True
        else:
            valid_start = False
        if nseq[coords[1] - 3:coords[1]] in STOP_CODONS:
            tx.att_tbl.pop('missing_stop_codon', None)
            valid_stop = True
        else:
            valid_stop = False
        if valid_start and valid_stop:
            tx.att_tbl['valid_ORF'] = 'True'
        tx.att_tbl['manual'] = 'True'
        return True, (tx.fid, coords[0], coords[1], cds_coords[0], cds_coords[1])
    return False, None

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
    
    wd = args.window if args.window != -1 else 300000000000

    if wd % 3 != 0:
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
    alt_start_stop = []
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        is_alt_found, coords = find_alt_start_stop(tx, cds_coords_tbl[tid], nfa[tid].seq, wd)
        if is_alt_found:
            alt_start_stop.append(coords) # tid, new_start, new_end, old_start, old_end
            cds_coords_tbl[tid] = (coords[1], coords[2])
    write_tup_lst(alt_start_stop, os.path.join(args.out_dir, 'alt_start_stop.csv'))
    print(tmessage(f"alternative start and/or stop found for {len(alt_start_stop)} transcripts", Mtype.PROG))

    # strategies to improve pident / resemblance with reference prots
    # 1. check containment
    # 2. check alt ORFs (optional)
    print(tmessage("checking reference orf containments", Mtype.PROG))
    contained = []
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        if tx.att_tbl['origin_ID'] in rslippage_tids: 
            tx.att_tbl['exception'] = 'ribosomal slippage'
            continue
        is_contained, coords = check_containment(tx, nfa[tid].seq, cds_coords_tbl[tid], \
                                            rfa[tx.att_tbl['origin_ID']].seq)
        if is_contained:
            contained.append(coords)
            cds_coords_tbl[tid] = (coords[1], coords[2])
    write_tup_lst(contained, os.path.join(args.out_dir, 'contained.csv'))
    print(tmessage(f"{len(contained)} containments detected and addressed", Mtype.PROG))

    if args.find_alt_orfs:
        alt_orfs = []
        for tid in tqdm(gan.txes):
            tx = gan.txes[tid]
            if tx.att_tbl['origin_ID'] in rslippage_tids: continue
            ao_found, coords = find_alt_orf(tid, nfa[tid].seq, rfa[tx.att_tbl['origin_ID']].seq, \
                                        cds_coords_tbl[tid], 0.9, 10)
            if ao_found:
                alt_orfs.append(coords)
                cds_coords_tbl[tid] = (coords[1], coords[2])
        write_tup_lst(alt_orfs, os.path.join(args.out_dir, 'alt_orfs.csv'))
        print(tmessage(f"{len(alt_orfs)} alternative orfs deteced and addressed", Mtype.PROG))
    
    print(tmessage("updating ref protein match and valid_ORFs tags", Mtype.PROG))
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        if tx.att_tbl['origin_ID'] in rslippage_tids: continue
        nseq = nfa[tid].seq
        cds_coords = cds_coords_tbl[tid]
        tsl_cseq = Seq.translate(nseq[cds_coords[0]:cds_coords[1]])
        if tsl_cseq[-1] == '*': tsl_cseq = tsl_cseq[:-1]
        rseq = rfa[tx.att_tbl['origin_ID']].seq
        if tsl_cseq == rseq:
            tx.att_tbl['matches_ref_protein'] = 'True'
    
    for gid in tqdm(gan.genes):
        gene = gan.genes[gid]
        vo_ctr = 0
        for tid in gene.children:
            tx = gan.txes[tid]
            if tx.att_tbl['valid_ORF'] == 'True':
                vo_ctr += 1
        gene.att_tbl['valid_ORFs'] = vo_ctr
    
    print(tmessage(f"saving results", Mtype.PROG))
    s = gan.to_str(args.fmt)

    with open(os.path.join(args.out_dir, f'swept.{args.fmt}'), 'w') as fh: fh.write(s)
    save_pth = os.path.join(args.out_dir, 'swept.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)