from tidy.utils import *
import pickle

def is_valid_stop(s) -> bool:
    last_codon = s[-3:].upper()
    if last_codon in STOP_CODONS:
        return True
    return False

def is_valid_start(s) -> bool:
    start_codon = s[:3].upper()
    if start_codon in START_CODONS:
        return True
    return False

def load_tsle(fn) -> dict:
    df = pd.read_csv(fn)
    tsle = dict()
    for _, row in df.iterrows():
        tsle[row['transcript_id']] = row['exception']
    return tsle

def check_start_and_stop(tx, cseq, tsle):
    val_start = is_valid_start(cseq)
    val_stop = is_valid_stop(cseq)
    if val_start and val_stop:
        tx.att_tbl['valid_ORF'] = 'True'
        tx.att_tbl.pop('missing_start_codon', None)
        tx.att_tbl.pop('missing_stop_codon', None)
        return
    else:
        ori_id = tx.att_tbl['origin_ID']
        if ori_id in tsle:
            if not val_start and tsle[ori_id] in START_TSLE:
                tx.att_tbl['exception'] = tsle[ori_id]
                val_start = True
            if not val_stop and tsle[ori_id] in STOP_TSLE: 
                tx.att_tbl['exception'] = tsle[ori_id]
                val_stop = True
            if val_start and val_stop:
                tx.att_tbl['valid_ORF'] = 'True'
                tx.att_tbl.pop('missing_start_codon', None)
                tx.att_tbl.pop('missing_stop_codon', None)
                return
        if not val_start:
            tx.att_tbl['missing_start_codon'] = 'True'
        if not val_stop:
            tx.att_tbl['missing_stop_codon'] = 'True'
        tx.att_tbl['valid_ORF'] = 'False'

# TODO: check if this is correct
def check_ptc(tx, pseq, tsle) -> bool:
    ori_id = tx.att_tbl['origin_ID']
    if '.' in pseq or '*' in pseq:
        if ori_id in tsle:
            if tsle[ori_id] in PTC_TSLE:
                tx.att_tbl['exception'] = tsle[ori_id]
                tx.att_tbl['valid_ORF'] = 'True'; return False
            else:
                tx.att_tbl['valid_ORF'] = 'False'; return True
        tx.att_tbl['valid_ORF'] = 'False'; return True
    return False


def main(args):
    if not is_pickled(args.in_file):
        print(tmessage("input file must a pickle object", Mtype.ERR))
        sys.exit(-1)
    
    print(tmessage("loading pickled gan", Mtype.PROG))
    with open(args.in_file, 'rb') as f: gan = pickle.load(f)

    if args.tsl_except:
        print(tmessage("loading translational exceptions", Mtype.PROG))
        if not os.path.exists(args.tsl_except):
            print(tmessage("tsl exception file not found", Mtype.ERR))
            sys.exit(-1)
        tsle = load_tsle(args.tsl_except)
    else:
        tsle = None
    
    if args.mane_file:
        print(tmessage("loading MANE transcript ids", Mtype.PROG))
        if not os.path.exists(args.mane_file):
            print(tmessage("MANE file not found", Mtype.ERR))
            sys.exit(-1)
        mane_tids = load_lst(args.mane_file)
    else:
        mane_tids = None

    cfa = pyfastx.Fasta(args.cds_file)
    pfa = pyfastx.Fasta(args.prot_file)
    rfa = pyfastx.Fasta(args.ref_file)

    print(tmessage("checking ORF validity, MANE status, and reference protein match states", Mtype.PROG))
    for tid in tqdm(gan.txes):
        tx = gan.txes[tid]
        check_start_and_stop(tx, cfa[tid].seq, tsle)
        has_ptc = check_ptc(tx, pfa[tid].seq, tsle)
        if has_ptc:
            tx.att_tbl['inframe_stop_codon'] = 'True'
        else:
            tx.att_tbl.pop('inframe_stop_codon', None)

        if tx.att_tbl['origin_ID'] in mane_tids:
            tx.att_tbl['tag'] = "MANE Select"
        if pfa[tid].seq == rfa[tx.att_tbl['origin_ID']].seq:
            tx.att_tbl['matches_ref_protein'] = 'True'
        else:
            tx.att_tbl['matches_ref_protein'] = 'False'
    
    print(tmessage("checking the number of valid ORFs", Mtype.PROG))
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

    with open(os.path.join(args.out_dir, f'labeled.{args.fmt}'), 'w') as fh: fh.write(s)
    save_pth = os.path.join(args.out_dir, 'labeled.pkl')
    with open(save_pth, 'wb') as f: pickle.dump(gan, f)