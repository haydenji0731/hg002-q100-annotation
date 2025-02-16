#!/usr/bin/env python

from utils import *
import name

def parse():
    parser = argparse.ArgumentParser(description="")
    subparsers = parser.add_subparsers(dest='module', help="")

    # assumes a sorted input gff file
    parser_name = subparsers.add_parser('name', help="")
    parser_name.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_name.add_argument('-o', '--out-file', type=str, help="", required=True)
    parser_name.add_argument('-p', '--prefix', type=str, help="", required=True)
    parser_name.add_argument('--ext', type=str, help="", required=False, default="gff")
    parser_name.add_argument('-clim', '--char-limit', type=int, help="", required=False, default=9)
    
    parser_init = subparsers.add_parser('init', help="")
    parser_init.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_init.add_argument('-o', '--out-dir', type=str, help="", required=True)
    args = parser.parse_args()
    return args

def main():
    args = parse()
    if args.module == 'name':
        print(tmessage(f'renaming features', Mtype.PROG))
        name.main(args)
    print(tmessage(f'finished', Mtype.PROG))

if __name__ == "__main__":
    main()