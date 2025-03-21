#!/usr/bin/env python

from tidy.utils import *
from tidy import name, build, trim, label, sweep

def parse():
    parser = argparse.ArgumentParser(description="")

    # TODO: add -i as a common arg
    parser.add_argument('--fmt', type=str, help="", required=False, \
                    default="gff", choices=['gff', 'gtf'])
    parser.add_argument('-o', '--out-dir', type=str, help="", required=True)
    parser.add_argument('-i', '--in-file', type=str, help="", required=True)

    subparsers = parser.add_subparsers(dest='module', help="")

    # assumes a sorted input gff file
    parser_name = subparsers.add_parser('name', help="")
    parser_name.add_argument('-fp', '--file-prefix', type=str, help="", \
                            default="renamed", required=False)
    parser_name.add_argument('-ip', '--id-prefix', type=str, help="", required=True)
    parser_name.add_argument('-clim', '--char-limit', type=int, help="", required=False, default=-1)
    

    parser_build = subparsers.add_parser('build', help="")
    parser_build.add_argument('--btype', action='store_true', default=False, required=False, help="")
    parser_build.add_argument('--schema', type=lambda s: s.split(','), required=False, \
                default=['ID', 'Parent', 'transcript_biotype'])
    
    parser_trim = subparsers.add_parser('trim', help="")
    parser_trim.add_argument('-c', '--cds-file', type=str, help="", required=True)
    parser_trim.add_argument('-r', '--ref-file', type=str, help="", required=True)
    parser_trim.add_argument('-d', '--min-diff', type=float, help="", required=False, default=5.0)

    parser_label = subparsers.add_parser('label', help="")
    parser_label.add_argument('-c', '--cds-file', type=str, help="", required=True)
    parser_label.add_argument('-r', '--ref-file', type=str, help="", required=True)
    parser_label.add_argument('-p', '--prot-file', type=str, help="", required=True)
    parser_label.add_argument('-e', '--tsl-except', type=str, help="", required=False)
    parser_label.add_argument('-m', '--mane-file', type=str, help="mane select", required=False)

    parser_sweep = subparsers.add_parser('sweep', help="")
    parser_sweep.add_argument('-w', '--window', type=int, help="", required=False, default=-1)
    parser_sweep.add_argument('-n', '--nucl-file', type=str, help="", required=True)
    parser_sweep.add_argument('-r', '--ref-file', type=str, help="", required=True)
    parser_sweep.add_argument('-e', '--slip-except', type=str, help="", required=False)
    parser_sweep.add_argument('--find-alt-orfs', action='store_true', default=False, required=False, help="")
    
    args = parser.parse_args()
    return args

def main():
    args = parse()
    check_dir(args.out_dir)
    # TODO: add an 'all' module that runs all of these ops
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
    elif args.module == 'label':
        print(tmessage(f'adding auxiliary information', Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "label_params.json")
        store_params(args, param_fn)
        label.main(args)
    elif args.module == 'sweep':
        print(tmessage(f'looking for CDS corrections', Mtype.PROG))
        param_fn = os.path.join(args.out_dir, "sweep_params.json")
        store_params(args, param_fn)
        sweep.main(args)
    print(tmessage(f'finished', Mtype.PROG))

if __name__ == "__main__":
    main()