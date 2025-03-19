from utils import *

def extend_upstream(s, cds_start, w) -> int:
    wd_end = max(0, cds_start - w) # respect tx ends
    for i in range(cds_start, wd_end, -3):
        if s[i - 3:i] in START_CODONS:
            return i - 3
    return -1

# TODO: finish implementation
def extend_downstream(s, cds_end, w) -> int:
    return -1

# TODO: finish implementation
def find_alt_start_stop(tx, cds_coords, nseq, w):
    if tx.att_tbl['valid_ORF'] == 'True': return None
    else:
        if 'inframe_stop_codon' in tx.att_tbl and \
            tx.att_tbl['inframe_stop_codon'] == 'True': return None
        if 'missing_start_codon' in tx.att_tbl and \
            tx.att_tbl['missing_start_codon'] == 'True':
            start_pos = extend_upstream(nseq, cds_coords[0], w)
            stop_pos = extend_downstream(nseq, cds_coords[0], w)
    return start_pos, stop_pos

# TODO: implement
def is_orf_contained():
    return

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

    nfa = pyfastx.Fasta(args.nucl_file)
    cds_coords_tbl = load_cds_coords(args.nucl_file)
    # assert len(cds_coords_tbl) == len(gan.txes) # sanity check

    for tid in gan.txes:
        tx = gan.txes[tid]
        find_alt_start_stop(tx, cds_coords_tbl[tid], nfa[tid].seq, args.window)
        

