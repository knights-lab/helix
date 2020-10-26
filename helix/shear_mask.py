#! /usr/bin/env python
import argparse
import datetime
import os
import csv
from collections import defaultdict
import numpy as np


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for blast mask',
        usage='strip_taxamap -f <input_fasta> -b <input_blast> -o <output_file>'
    )
    parser.add_argument('-f', '--fasta', help='Set the directory path of the input fasta file', required=True)
    parser.add_argument('-b', '--blast', help="Set the path of the input blast file", required=True)
    parser.add_argument('-r', '--read_length', help='Set the read length of output reads', required=True, type=int)
    parser.add_argument('-o', '--output', help='Set the directory path of the output masked fasta file (default: cwd)', required=True)
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


def gen_blast_reads_hits(inf_path: str, shear_length: int):
    with open(inf_path) as inf:
        csv_reader = csv.reader(inf, delimiter="\t")
        for row in csv_reader:
            reference = "_".join(row[0].split("_")[:-1])
            hit_begin = int(row[0].split("_")[-1])
            hit_end = hit_begin + shear_length
            yield reference, (hit_begin, hit_end)


def build_mask_dict(inf_path: str, shear_length: int) -> dict:
    dd = defaultdict(list)
    gen = gen_blast_reads_hits(inf_path, shear_length)
    for row in gen:
        dd[row[0]].append(row[1])
    return dd


def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = ''
    for line in fh:
        if line[0] == ">":
            if title:
                yield title.strip(), data
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield title.strip(), data


def mask_fasta(gen_fasta, mask_dict):
    for title, data in gen_fasta:
        np_data = np.fromiter(data, dtype='S1', count=len(data))
        for mask in mask_dict[title]:
            np_data[mask[0] - 1:mask[1] - 1] = "N"
        yield title, np_data.tobytes().decode()


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    dd_mask_hits = build_mask_dict(args.blast, args.read_length)

    with open(args.fasta, "rb") as inf:
        gen_fasta = read_fasta(inf)
        gen_mask_fasta = mask_fasta(gen_fasta, dd_mask_hits)
        with open(args.output, "w") as outf:
            for header, seq in gen_mask_fasta:
                outf.write(f">{header}\n{seq}\n")

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
