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
    parser.add_argument('-r', '--read_length', help='Set the read length of output reads', required=True)
    parser.add_argument('-s', '--step_size', help='Set the step size of the output file', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output fasta (default: cwd)', default=os.path.join(os.getcwd(), "shear.fna"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def shear_fasta(inf: str, outf: str, read_length: int, step_size: int):
    with open(outf, "w") as outf:
        for line in open(inf):
            if line.startswith('>'):
                header = line.strip()
                header = header.replace(".", "_", 1)
                header_words = header.split()
                header_id = header_words[0]
                header_comments = ''
                if len(header_words) > 1:
                    header_comments = ' ' + ' '.join(header_words[1:])
            else:
                seq = line.strip()
                n = len(seq)
                start_pos = 0
                while start_pos < n:
                    outf.write(header_id + "_%09d" % (start_pos) + header_comments + "\n")
                    if start_pos + read_length > n:
                        start_pos = n - read_length
                    outf.write(seq[start_pos:(start_pos + read_length)] + "\n")
                    if start_pos == (n - read_length):
                        start_pos = n
                    else:
                        start_pos += step_size


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.abspath(os.path.join(args.output))
    os.makedirs(outdir, exist_ok=True)

    shear_fasta(args.fasta, args.output, args.read_length, args.step_size)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
