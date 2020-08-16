#! /usr/bin/env python
import argparse
import datetime
import os
import csv


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for strip_taxamap',
        usage='strip_taxamap <input_blast> -i <input> -o <output_file>'
    )
    parser.add_argument('-i', '--input', help='Set the directory path of the input mapping file', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output taxamap file (default: cwd)', required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def strip_dangling_taxa(taxa_str: str) -> str:
    taxa_l = taxa_str.split(";")
    for ix in range(len(taxa_l) - 1, -1, -1):
        taxa_str = taxa_l[ix]
        if len(taxa_str) > 3 and not taxa_str.endswith("__"):
            break
    taxa_str = ";".join(taxa_l[:ix + 1])
    return taxa_str


def gen_strip_dangling_taxa(inf_path):
    with open(inf_path) as inf:
        csv_reader = csv.reader(inf, delimiter="\t")
        for row in csv_reader:
            taxa_str = row[1]
            strip_str = strip_dangling_taxa(taxa_str)
            strip_str = strip_str.replace(" ", "_")
            yield row[0], strip_str


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    gen = gen_strip_dangling_taxa(args.input)

    with open(args.output, "w") as outf:
        for row in gen:
            outf.write("\t".join(row) + "\n")

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
