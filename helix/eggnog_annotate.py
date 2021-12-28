import csv
import os.path
import pathlib
from collections import defaultdict
import argparse

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for filter_fasta',
    )
    parser.add_argument('-t', '--taxa_map', help='Set the path to the input taxamap', required=True)
    parser.add_argument('-a', '--annotation', help='Set the path to the csv annotation file', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.path.join(os.getcwd() + "results.txt"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def get_taxamap(inf_path: pathlib.Path) -> dict:
    taxa_map = dict()
    with open(inf_path) as inf:
        csv_reader = csv.reader(inf, delimiter="\t")
        for row in csv_reader:
            taxa_map[row[0]] = row[1]
    return taxa_map


def get_annotations(taxa_map: dict, annotations: pathlib.Path) -> dict:
    results = defaultdict(list)
    num_in = 0
    num_out = 0
    with open(annotations) as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            else:
                row = line.split("\t")
                query_list = row[0].split(".")
                query = "_".join(query_list[1:3]) + "." + query_list[3]
                if query in taxa_map:
                    num_in += 1
                    ko = row[11]
                    if ko.startswith("ko:"):
                        results[query].append(ko.split(":")[1])
                else:
                    num_out += 1
                    # raise Exception(f"{query} not in Map")
    return results

def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    taxa_map = get_taxamap(args.taxa_map)
    results = get_annotations(taxa_map, args.annotation)

    with open(args.output, "w") as outf:
        for k, v in results.items():
            outf.write(k + "\t" + "\t".join(v) + "\n")

if __name__ == '__main__':
    main()
