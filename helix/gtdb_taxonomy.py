#! /usr/bin/env python
# Shears fasta into substrings of length w step d
# usage
# shear_db.py db.fasta 100 50
import argparse
import datetime
import os
import sys


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
    headers_fasta = set()
    headers_map = set()
    # we need assert that all headers are unique in both the taxamap
    # make sure that all headers are found in the taxamap
    with open(fasta_inf) as inf:
        for line in inf:
            if line.startswith(">"):
                header = line.rstrip()[1:]
                if header in headers_fasta:
                    sys.stderr.write(f"WARNING: the header {header} found twice in the input FASTA {fasta_inf}")
                headers_fasta.add(header)
    with open(outf, "w") as outf:
        with open(tax_inf) as inf:
            for line in inf:
                line = line.rstrip()
                header, tax = line.split("\t")
                tax = tax.replace(" ", "_")
                tax = strip_dangling_taxa(tax)
                # header = "_".join(header.split("_")[1:])
                if header in headers_fasta:
                    outf.write(f"{header}\t{tax}\n")
                    if header in headers_map:
                       sys.stderr.write(f"WARNING: the header {header} found twice in the input map {tax_inf}")
                    headers_map.add(header)
    headers_difference = headers_fasta.difference(headers_map)
    if len(headers_difference) > 0:
        sys.stderr.write(f"WARNING: the headers {headers_difference} not found in the input map")


def strip_dangling_taxa(taxa_str: str) -> str:
    taxa_l = taxa_str.split(";")
    for ix in range(len(taxa_l) - 1, -1, -1):
        taxa_str = taxa_l[ix]
        if len(taxa_str) > 3 and not taxa_str.endswith("__"):
            break
    taxa_str = ";".join(taxa_l[:ix + 1])
    return taxa_str


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
