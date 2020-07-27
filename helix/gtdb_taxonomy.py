#! /usr/bin/env python
# Shears fasta into substrings of length w step d
# usage
# shear_db.py db.fasta 100 50
import argparse
import datetime
import os


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for shear_db',
        usage='shear_db %s -b <input_blast> -f <input_fastq> -o <output_file>'
    )
    parser.add_argument('-f', '--fasta', help='Set the directory path of the input fasta', required=True)
    parser.add_argument('-t', '--taxa_table', help='Set the directory path of the input taxafile', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output taxatable (default: cwd)', default=os.path.join(os.getcwd(), "shear.fna"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def gtdb_taxonomy(fasta_inf: str, outf: str, tax_inf: str):
    headers = set()
    with open(fasta_inf) as inf:
        for line in inf:
            if line.startswith(">"):
                headers.add(line.rstrip()[1:])
    with open(outf, "w") as outf:
        with open(tax_inf) as inf:
            for line in inf:
                line = line.rstrip()
                header, tax = line.split("\t")
                header = "_".join(header.split("_")[1:])
                if header in headers:
                    outf.write(f"{header}\t{tax}\n")


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    gtdb_taxonomy(args.fasta, args.output, args.taxa_table)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
