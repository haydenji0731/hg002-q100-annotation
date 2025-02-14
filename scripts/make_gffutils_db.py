#!/usr/bin/env python

from commons import *
import gffutils
import sys

def load_id_spec(fn):
    fh = open(fn, 'r')
    id_spec = dict()
    for ln in fh:
        temp = ln.strip().split(",")
        if len(temp) < 2:
            print(tmessage(f'id_spec should contain 2 columns', Mtype.ERR))
            sys.exit(-1)
        if temp[0] == '*':
            key_atts = []
            for i in range(1, len(temp), 1):
                key_atts.append(temp[i])
            return key_atts
        else:
            key_atts = []
            for i in range(len(temp)):
                if i == 0:
                    feature = temp[i]
                else:
                    key_atts.append(temp[i])
            if len(key_atts) == 1:
                id_spec[feature] = key_atts[0]
            else:
                id_spec[feature] = key_atts
    fh.close()
    return id_spec

def main(in_fn, out_fn, spec_fn):
    id_spec = load_id_spec(spec_fn)
    _ = gffutils.create_db(in_fn, dbfn=out_fn, force=True, \
                        disable_infer_genes=True, disable_infer_transcripts=True, \
                        id_spec=id_spec)
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])