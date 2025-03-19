from utils import *
from Bio.Seq import Seq
import psa
import copy

# def load_tsle(fn) -> dict:
#     df = pd.read_csv(fn)
#     tsle = dict()
#     for _, row in df.iterrows():
#         tsle[row['transcript_id']] = row['exception']
#     return tsle

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


def start_or_end(cseq, r_pseq, l) -> int:
    temp = cseq[l:]
    temp_pseq = Seq.translate(temp)
    aln = psa.water(moltype='prot', qseq=temp_pseq, sseq=r_pseq)
    psim_1 = aln.pidentity
    temp = cseq[:-l]
    temp_pseq = Seq.translate(temp)
    aln = psa.water(moltype='prot', qseq=temp_pseq, sseq=r_pseq)
    psim_2 = aln.pidentity
    if psim_1 > psim_2:
        return 0
    elif psim_2 > psim_1:
        return 1
    return -1 # not undetermined

def check_cds_integrity(cds_chain):
    is_valid = True
    for cds in cds_chain:
        if cds.end < cds.start:
            is_valid = False
            break
    return is_valid

def shift_cds(tx, l, side):
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
    if not args.in_file.lower().endswith('.pkl'):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: in_gan = pickle.load(f)

    print(tmessage("calculating CDS lengths", Mtype.PROG))
    clen_tbl = {tid: calc_clen(tx) for tid, tx in in_gan.txes.items()}
    for tx in in_gan.txes.values(): 
        if zero_start_frame(tx): clen_tbl[tx.fid] = calc_clen(tx)

    pfa = pyfastx.Fasta(args.prot_file)
    nfa = pyfastx.Fasta(args.nucl_file)
    cfa = pyfastx.Fasta(args.cds_file)
    r_pfa = pyfastx.Fasta(args.ref_file)
    plen_tbl = {x.name: len(x.seq) for x in pfa}

    manual = []
    shifted = []
    ctr = 0
    for tid in tqdm(clen_tbl):
        mod = clen_tbl[tid] % 3
        tx = in_gan.txes[tid]
        if mod != 0:
            ctr += 1
            res = start_or_end(cfa[tid].seq, r_pfa[tx.att_tbl['origin_ID']].seq, mod)
            if res == -1:
                manual.append(tid); continue
            else:
                shifted.append(tid)
            # else:
            #     shifted_cds_chain = shift_cds(tx, mod, res)
            #     tx.chains[1] = shifted_cds_chain
            #     tx.sort_chain(1)
    print(tmessage(f"{ctr} coding transcripts with CDS len indivisible by 3", Mtype.PROG))
    print(tmessage(f"{len(shifted)} / {ctr} trimmed", Mtype.PROG))
    print(tmessage(f"{len(manual)} / {ctr} needs manual inspection", Mtype.PROG))
    


    # if args.tsl_except:
    #     if not os.path.exists(args.tsl_except):
    #         print(tmessage("tsl exception file not found", Mtype.ERR))
    #         sys.exit(-1)
    #     tsle = load_tsle(args.tsl_except)
    
