import argparse
import os
import csv
import datetime
from itertools import zip_longest


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for filter_fastq',
        usage='filter_fastq %s -b <input_blast> -f <input_fastq> -o <output_folder>'
    )
    parser.add_argument('-f', '--fastq', help='Set the directory path of the fastq', required=True)
    parser.add_argument('-b', '--blast', help='Set the directory path of the blast alignment', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())

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


def read_fastq(fh):
    # Assume linear FASTQS
    count = 0
    for (title, sequence, garbage, qualities) in grouper(fh, 4, fillvalue=""):
        count += 1
        # Record begins
        if title[0] != '@':
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records - Title line does not begin with "@" symbol - error on line' + str(count) + '.')
        title = title[1:].strip()
        sequence = sequence.strip()
        count += 1
        garbage = garbage.strip()
        count += 1
        if garbage[0] != '+':
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records - strand line does not contain "+" symbol on line' + str(count) + '.')
        qualities = qualities.strip()
        count += 1
        if len(qualities) != len(sequence):
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records - Sequence length does not equal quality score length on line ' + str(count) + '.')
        yield title, sequence, qualities


def filter_reads(fp_fastq: str, outdir: str, filter_l: set) -> None:
    f_reads = 0
    k_reads = 0

    base = ".".join(os.path.basename(fp_fastq).split(".")[:-1])

    with open(fp_fastq) as inf_fastq:
        with open(os.path.join(outdir, f"{base}.fq"), "w") as outf_insample:
            with open(os.path.join(outdir, f"{base}.filtered.fq"), "w") as outf_outsample:
                fastq_gen = read_fastq(inf_fastq)
                for ix, (title, sequence, qualities) in enumerate(fastq_gen):
                    name = f"{base}_{ix}"
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

    filter_l = build_read_filter_set(args.blast)

    filter_reads(args.fastq, outdir, filter_l)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
