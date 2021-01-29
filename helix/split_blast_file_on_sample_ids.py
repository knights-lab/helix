#! /usr/bin/env python
import argparse
import datetime
import os
from itertools import zip_longest
from collections import defaultdict
import csv


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for split_sequence_file_on_sample_ids',
    )
    parser.add_argument('-i', '--input', help='Set the directory path of the input blast file', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output folder (default: cwd)', default=os.getcwd())
    parser.add_argument('-b', '--buffer', help="the buffer size", type=int, default=1_000_000)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def read_blast(inf):
    inf_csv = csv.reader(inf, delimiter="\t")
    for row in inf_csv:
        yield row[0], "\t".join(row).strip() + "\n"


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # taken from itertools recipes: https://docs.python.org/3/library/itertools.html#itertools-recipes
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def split_sequence_file(input, output_dir, buffer=1_000_000):
    with open(input) as inf:
        for group in grouper(read_blast(inf), buffer):
            d_group = defaultdict(str)
            for header, row in group:
                sample_id = header.split()[0].split("_")[0]
                d_group[sample_id] = row
            for k, v in d_group.items():
                outfile = os.path.join(output_dir, f"{k}.b6")
                with open(outfile, "a+") as outfile:
                    outfile.write(v)


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    # outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(args.output, exist_ok=True)

    split_sequence_file(args.input, args.output, buffer=args.buffer)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
