import argparse
import os
import csv
import datetime
from itertools import zip_longest


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for filter_fasta',
    )
    parser.add_argument('-f', '--fasta', help='Set the directory path of the fasta', required=True)
    parser.add_argument('-a', '--alignment', help='Set the directory path of the alignment tsv', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('--original', default=False, action='store_true', help='Set if the sequence names are the original sequence names.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def build_read_filter_set(fp_blast: str) -> set:
    filter_l = set()
    with open(fp_blast) as inf_blast:
        reader = csv.reader(inf_blast, delimiter="\t")
        for row in reader:
            # add title to filter
            filter_l.add(row[0])
    return filter_l


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # taken from itertools recipes: https://docs.python.org/3/library/itertools.html#itertools-recipes
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = ""
    data = ""
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


def filter_reads(fp_fasta: str, outdir: str, filter_l: set, change_flag: bool = False) -> None:
    f_reads = 0
    k_reads = 0

    base = ".".join(os.path.basename(fp_fasta).split(".")[:-1])

    with open(fp_fasta) as inf_fasta:
        with open(os.path.join(outdir, f"{base}.fq"), "w") as outf_insample:
            with open(os.path.join(outdir, f"{base}.filtered.fq"), "w") as outf_outsample:
                fasta_gen = read_fasta(inf_fasta)
                for ix, (title, sequence, qualities) in enumerate(fasta_gen):
                    if not change_flag:
                        name = f"{base}_{ix}"
                    else:
                        name = title.split()[0]
                    if name in filter_l:
                        outf_outsample.write(f"@{title}\n{sequence}\n+\n{qualities}\n")
                        f_reads += 1
                    else:
                        outf_insample.write(f"@{title}\n{sequence}\n+\n{qualities}\n")
                        k_reads += 1
    print(f"Filtered reads:\t{f_reads}\nKept reads:\t{k_reads}")


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.abspath(os.path.join(args.output))
    os.makedirs(outdir, exist_ok=True)

    filter_l = build_read_filter_set(args.alignment)

    filter_reads(args.fasta, outdir, filter_l, change_flag=args.original)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
