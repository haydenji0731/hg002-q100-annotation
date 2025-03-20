from utils import *
import psa

def calc_clen(x) -> int:
    clen = 0
    if not x.is_chain_init(): return -1
    for c in x.chains[1]:
        clen += c.end - c.start + 1
    return clen

def zero_start_frame(x) -> bool:
    assert x.sorted[1] # sanity check
    start_frame = x.chains[1][0].frame
    assert start_frame is not None
    if start_frame == 0:
        return False
    if x.strand == '+':
        x.chains[1][0].start += start_frame
    else:
        x.chains[1][0].end -= start_frame
    x.chains[1][0].frame = 0 # zero out the starting frame
    return True


def start_or_end(cseq, r_pseq, l):
    temp = cseq[l:]
    temp_pseq = Seq.translate(temp)
    aln = psa.water(moltype='prot', qseq=temp_pseq, sseq=r_pseq)
    nident_1 = aln.nidentity
    temp = cseq[:-l]
    temp_pseq = Seq.translate(temp)
    aln = psa.water(moltype='prot', qseq=temp_pseq, sseq=r_pseq)
    nident_2 = aln.nidentity
    if nident_1 > nident_2:
        return 0, nident_1, nident_2
    elif nident_2 > nident_1:
        return 1, nident_1, nident_2
    return -1, nident_1, nident_2

def check_cds_integrity(cds_chain):
    is_valid = True
    for cds in cds_chain:
        if cds.end < cds.start:
            is_valid = False
            break
    return is_valid

def shorten_cds(tx, l, side):
    new_cds_chain = copy.deepcopy(tx.chains[1])
    remaining = l
    if tx.strand == '+':
        if side == 0: # start
            while remaining > 0:
                curr_clen = new_cds_chain[0].end - new_cds_chain[0].start + 1
                if remaining >= curr_clen:
                    new_cds_chain = new_cds_chain[1:]
                    remaining -= curr_clen
                else:
                    new_cds_chain[0].start += remaining
                    remaining -= remaining # sets to zero
        else:
            while remaining > 0:
                curr_clen = new_cds_chain[-1].end - new_cds_chain[-1].start + 1
                if remaining >= curr_clen:
                    new_cds_chain = new_cds_chain[:-1]
                    remaining -= curr_clen
                else:
                    new_cds_chain[-1].end -= remaining
                    remaining -= remaining
    else:
        if side == 0: # start
            while remaining > 0:
                curr_clen = new_cds_chain[0].end - new_cds_chain[0].start + 1
                if remaining >= curr_clen:
                    new_cds_chain = new_cds_chain[1:]
                    remaining -= curr_clen
                else:
                    new_cds_chain[0].end -= remaining
                    remaining -= remaining
        else:
            while remaining > 0:
                curr_clen = new_cds_chain[-1].end - new_cds_chain[-1].start + 1
                if remaining >= curr_clen:
                    new_cds_chain = new_cds_chain[:-1]
                    remaining -= curr_clen
                else:
                    new_cds_chain[-1].start += remaining
                    remaining -= remaining
    # sanity check
    assert check_cds_integrity(new_cds_chain)
    return new_cds_chain

def main(args):
    if not is_pickled(args.in_file):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: gan = pickle.load(f)

    print(tmessage("calculating CDS lengths", Mtype.PROG))
    clen_tbl = {tid: calc_clen(tx) for tid, tx in gan.txes.items()}
    for tx in gan.txes.values(): 
        if zero_start_frame(tx): clen_tbl[tx.fid] = calc_clen(tx)

    cfa = pyfastx.Fasta(args.cds_file)
    r_pfa = pyfastx.Fasta(args.ref_file)

    manual = []
    start_shifted = []
    aln_shifted = []
    ctr = 0
    for tid in tqdm(clen_tbl):
        mod = clen_tbl[tid] % 3
        tx = gan.txes[tid]
        if mod != 0:
            ctr += 1
            cseq = cfa[tid].seq
            if cseq[:3] in START_CODONS:
                start_shifted.append((tid, str(mod)))
                shifted_cds_chain = shorten_cds(tx, mod, 1) # trim at the end
                tx.chains[1] = shifted_cds_chain
                tx.sort_chain(1)
                tx.att_tbl['manual'] = 'True'
            else:
                side, nident_start, nident_end = start_or_end(cseq, r_pfa[tx.att_tbl['origin_ID']].seq, mod)
                if side == -1:
                    manual.append((tid, str(mod), str(nident_start), str(nident_end)))
                    tx.att_tbl['manual'] = 'False'
                else:
                    aln_shifted.append((tid, str(mod), str(side), str(nident_start), str(nident_end)))
                    shifted_cds_chain = shorten_cds(tx, mod, side)
                    tx.chains[1] = shifted_cds_chain
                    tx.sort_chain(1)
                    tx.att_tbl['manual'] = 'True'
        else:
            tx.att_tbl['manual'] = 'False'

    print(tmessage(f"{ctr} coding transcripts with CDS len indivisible by 3", Mtype.PROG))
    print(tmessage(f"{len(start_shifted) + len(aln_shifted)} / {ctr} trimmed", Mtype.PROG))
    print(tmessage(f"{len(manual)} / {ctr} needs manual inspection", Mtype.PROG))


    print(tmessage(f"saving results", Mtype.PROG))
    write_tup_lst(start_shifted, os.path.join(args.out_dir, 'start_shifted.csv'))
    write_tup_lst(aln_shifted, os.path.join(args.out_dir, 'aln_shifted.csv'))
    write_tup_lst(manual, os.path.join(args.out_dir, 'manual_shifts.csv'))
    s = gan.to_str(args.fmt)

    with open(os.path.join(args.out_dir, f'trimmed.{args.fmt}'), 'w') as fh: fh.write(s)
    save_pth = os.path.join(args.out_dir, 'trimmed.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)