from datetime import datetime
from enum import Enum
import sys
import pandas as pd
import argparse
import pickle
import os
import json
import pyfastx
from tqdm import tqdm
import copy
from Bio.Seq import Seq
from subprocess import call

RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
RESET = '\033[0m'

# TODO: confirm these tsl exceptions
START_TSLE = ['aa:Met', 'aa:Leu']
STOP_TSLE = ['aa:TERM']
PTC_TSLE = ['aa:Sec']

START_CODONS = ['ATG', 'CTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

class Mtype(Enum):
    PROG = (GREEN, "PROGRESS")
    ERR = (RED, "ERROR")
    WARN = (YELLOW, "WARNING")

def tmessage(s, mtype) -> str:
    if mtype not in Mtype:
        raise Exception("Error while printing message")
    return f"{datetime.now()} {mtype.value[0]}{mtype.value[1]}{RESET} {s}"

def split_ln(ln, sep):
    clean_ln = ln.strip()
    fields = clean_ln.split(sep)
    clean_fields = [it.strip() for it in fields]
    return clean_fields

def to_dot(v):
    out = v
    if v is None:
        out = '.'
    return out

def copy_o(in_o, oid, format, feature):
    out_o = gLine(None, format)
    out_o.ctg = in_o.ctg
    out_o.src = in_o.src
    out_o.start = in_o.start
    out_o.end = in_o.end
    out_o.strand = in_o.strand
    out_o.frame = in_o.frame
    out_o.score = in_o.score
    out_o.feature = feature
    out_o.attributes["ID"] = oid
    return out_o

def write_lst(lst, fn):
    with open(fn, 'w') as fh:
        for x in lst:
            fh.write(f'{x}\n')

def load_lst(fn) -> list:
    lst = []
    with open(fn, 'r') as fh:
        for ln in fh:
            x = ln.strip()
            lst.append(x)
    return lst

def write_tup_lst(lst, fn):
    with open(fn, 'w') as fh:
        for x in lst:
            fh.write(f'{",".join(map(str, x))}\n')

def check_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def store_params(args, fn):
    with open(fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

def is_pickled(fn):
    return fn.lower().endswith('.pkl')

# custom classes
class gLine():
    def __init__(self, ln, format):
        if ln is None:
            self.init_empty()
        else:
            if format.lower() == 'gtf':
                kv_sep = ' '
                quoted = True
            elif format.lower() == 'gff':
                kv_sep = '='
                quoted = True
            else:
                print(tmessage("invalid file format", Mtype.ERR))
                sys.exit(-1)
            fields = split_ln(ln, sep='\t')
            if len(fields) != 9:
                print(tmessage("a line must have 9 columns", Mtype.ERR))
                sys.exit(-1)

            self.ctg = fields[0]
            self.src = fields[1]
            self.feature = fields[2]
            self.start = int(fields[3])
            self.end = int(fields[4])
            self.score = float(fields[5]) if fields[5] != '.' else None
            self.strand = fields[6]
            self.frame = fields[7] if fields[7] != '.' else None

            temp = split_ln(fields[8], sep=';')
            self.attributes = dict()
            for att in temp:
                kv_pair = att.split(kv_sep)
                if len(kv_pair) < 2:
                    continue
                k = kv_pair[0]
                v = kv_pair[1]
                if quoted:  v = v.replace('"', '')
                self.attributes[k] = v
    
    def init_empty(self) -> None:
        self.ctg = None
        self.src = None
        self.feature = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.frame = None
        self.attributes = dict()
    
    def to_gStr(self, format) -> str:
        gStr = f'{self.ctg}\t{self.src}\t{self.feature}\t{self.start}\t{self.end}'
        gStr += f'\t{to_dot(self.score)}\t{self.strand}\t{to_dot(self.frame)}\t'
        temp = list(self.attributes.items())
        kv_sep = "=" if format.lower() == "gff" else " "
        temp = [f'{x[0]}{kv_sep}{x[1]}' for x in temp]
        gStr += ';'.join(temp) + '\n'
        return gStr

class gFeat():
    def __init__(self, row: pd.Series, att_sep: str, schema: list):
        self.ctg = row['ctg']
        self.src = row['src']
        self.type = row['type']
        self.start = int(row['start'])
        self.end = int(row['end'])
        self.score = int(row['score']) if row['score'] != '.' else None
        self.strand = row['strand']
        self.frame = int(row['frame']) if row['frame'] != '.' else None
        self.att_tbl = self.att2dict(row['attributes'], att_sep)
        self.children = set()

        # need flexibility
        self.fid = self.att_tbl[schema[0]] if schema[0] in self.att_tbl else None
        self.parent = self.att_tbl[schema[1]] if schema[1] in self.att_tbl else None
        self.ftype = self.att_tbl[schema[2]] if schema[2] in self.att_tbl else None
        self.chains = []
        self.sorted = []
    
    def att2dict(self, s, sep) -> dict:
        fields = s.strip().split(';')
        att = dict()
        for x in fields:
            temp = x.strip().split(sep)
            if len(temp) < 2: continue
            k = temp[0]
            v = temp[1]
            att[k] = v
        return att

    # tx-specific fns
    def is_chain_init(self) -> bool:
        return len(self.chains) > 0
    
    def init_chain(self) -> bool:
        self.chains.append([])
        self.chains.append([])
        self.sorted.append(False)
        self.sorted.append(False)
        assert len(self.chains) == 2
    
    def print_chain(self, i):
        coords = [(x.start, x.end) for x in self.chains[i]]
        print(' '.join([f'({x}, {y})' for x, y in coords]))
    
    def append2chain(self, i, el):
        self.chains[i].append(el)
        self.sorted[i] = False
    
    def sort_chain(self, i):
        if self.strand == '+':
            self.chains[i] = sorted(self.chains[i], \
                                key=lambda x: x.start, reverse=False)
        else:
            self.chains[i] = sorted(self.chains[i], \
                                key=lambda x: x.end, reverse=True)
        self.sorted[i] = True
    
    def set_chain(self, i, chain):
        self.chains[i] = chain
        self.sorted[i] = False
    
    def get_chain_lens(self):
        return len(self.chains[0]), len(self.chains[1])
    
    # gene-specific fns
    def add_child(self, el):
        self.children.add(el)
    
    def has_child(self, el) -> bool:
        return el in self.children
    
    def __eq__(self, o) -> bool:
        if self.start != o.start or self.end != o.end: return False
        if self.strand != o.strand: return False
        if self.is_chain_init() and o.is_chain_init():
            for i in range(2):
                if not o.sorted[i] or not self.sorted[i]: print("can't perform eq on unsorted chains")
                if len(self.chains[i]) != len(o.chains[i]): return False
                for j in range(len(self.chains[i])):
                    if self.chains[0][j].start != o.chains[0][j].start or \
                        self.chains[0][j].end != o.chains[0][j].end: return False
        return True
    
    def to_str(self, sep) -> str:
        s = f'{self.ctg}\tLiftoff\t{self.type}\t{self.start}\t{self.end}\t'
        frame = '.' if self.frame is None else self.frame
        s += f'.\t{self.strand}\t{frame}\t'
        if 'ID' in self.att_tbl:
            s += f'ID{sep}{self.att_tbl["ID"]};'
        for k in self.att_tbl:
            if k != 'ID':
                s += f'{k}{sep}{self.att_tbl[k]};'
        return s
    
class gAn():
    def __init__(self, gene_l: list, tx_l: list):
        self.genes = {x.fid: x for x in gene_l}
        self.txes = {x.fid: x for x in tx_l}
    
    def to_str(self, fmt) -> str:
        sep = "=" if fmt.lower() == 'gff' else ' '
        gene_ctr = 0
        tx_ctr = 0
        out_lns = []
        for gid in self.genes:
            gene = self.genes[gid]
            out_lns.append(gene.to_str(sep))
            gene_ctr += 1
            for tid in gene.children:
                tx = self.txes[tid]
                out_lns.append(tx.to_str(sep))
                for exon in tx.chains[0]:
                    out_lns.append(exon.to_str(sep))
                for cds in tx.chains[1]:
                    out_lns.append(cds.to_str(sep))
                tx_ctr += 1
        assert len(self.genes) == gene_ctr
        assert len(self.txes) == tx_ctr
        s = '\n'.join(out_lns) + '\n'
        return s