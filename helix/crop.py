#! /usr/bin/env python
# Cuts sequencing lengths to min and max
import argparse
import datetime
import os

from helix.utils.utils import read_fasta


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for shear_db',
        usage='shear_db %s -b <input_blast> -f <input_fastq> -o <output_file>'
    )
    parser.add_argument('-f', '--fasta', help='Set the directory path of the input fasta', required=True)
    parser.add_argument('--max_length', help='Set the max length of the sequencing reads', required=True, type=int)
    parser.add_argument('--min_length', help='Set the min length of the sequencing reads', required=True, type=int)
    parser.add_argument('-o', '--output', help='Set the directory path of the output taxatable (default: cwd)', default=os.path.join(os.getcwd(), "shear.fna"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def crop(fasta_in, fasta_out, min_length: int, max_length: int):
    with open(fasta_in) as inf:
        with open(fasta_out, "w") as outf:
            for header, seq in read_fasta(inf):
                if len(seq) >= min_length:
                    outf.write(f">{header}\n{seq[:max_length]}\n")


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    if args.max_length < args.min_length:
        raise Exception(f"Max length {args.max_length} must be greater than min length {args.min_length}")

    crop(args.fasta, args.output, args.min_length, args.max_length)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
