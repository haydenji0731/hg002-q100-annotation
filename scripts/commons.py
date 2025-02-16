from datetime import datetime
from enum import Enum
import sys
import os

RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
RESET = '\033[0m'

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