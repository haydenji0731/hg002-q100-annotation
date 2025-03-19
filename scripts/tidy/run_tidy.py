#!/usr/bin/env python

from utils import *
import name, build, trim

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--ext', type=str, help="", required=False, default="gff")
    parser.add_argument('-o', '--out-dir', type=str, help="", required=True)
    subparsers = parser.add_subparsers(dest='module', help="")

    # assumes a sorted input gff file
    parser_name = subparsers.add_parser('name', help="")
    parser_name.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_name.add_argument('-fp', '--file-prefix', type=str, help="", \
                            default="renamed", required=False)
    parser_name.add_argument('-ip', '--id-prefix', type=str, help="", required=True)
    parser_name.add_argument('-clim', '--char-limit', type=int, help="", required=False, default=-1)
    

    parser_build = subparsers.add_parser('build', help="")
    parser_build.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_build.add_argument('--btype', action='store_true', default=False, required=False, help="")
    parser_build.add_argument('--schema', type=lambda s: s.split(','), required=False, \
                default=['ID', 'Parent', 'transcript_biotype'])
    
    parser_trim = subparsers.add_parser('trim', help="")
    parser_trim.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_trim.add_argument('-p', '--prot-file', type=str, help="", required=True)
    parser_trim.add_argument('-n', '--nucl-file', type=str, help="", required=True)
    parser_trim.add_argument('-c', '--cds-file', type=str, help="", required=True)
    parser_trim.add_argument('-r', '--ref-file', type=str, help="", required=True)
    parser_trim.add_argument('-e', '--tsl-except', type=str, help="", required=False)
    args = parser.parse_args()
    return args

def main():
    args = parse()
    check_dir(args.out_dir)
    if args.module == 'name':
        print(tmessage(f'renaming features', Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "name_params.json")
        store_params(args, param_fn)
        name.main(args)
    elif args.module == 'build':
        print(tmessage(f'building gan database', Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "build_params.json")
        store_params(args, param_fn)
        build.main(args)
    elif args.module == 'trim':
        print(tmessage(f'trimming CDS ends', Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "trim_params.json")
        store_params(args, param_fn)
        trim.main(args)
    print(tmessage(f'finished', Mtype.PROG))

if __name__ == "__main__":
    main()