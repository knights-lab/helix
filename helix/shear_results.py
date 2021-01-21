#! /usr/bin/env python
# usage: python me.py \
#  alignment.burst.otu.txt db.tax sheared_bayes.txt
import csv
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import argparse
import datetime
import os


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for shear_db',
        usage='shear_db %s -b <input_blast> -f <input_fastq> -o <output_file>'
    )
    parser.add_argument('-t', '--taxonomy_table', help='Set the directory path of the input taxonomy table', required=True)
    parser.add_argument('-m', '--taxonomy_map', help='Set the directory path of the input taxonomy map', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output sheared_bayes file (default: cwd)', default=os.path.join(os.getcwd(), "sheared_bayes.txt"))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ')
    return parser


def index_lca(str1: str, str2: str) -> int:
    for i, (s1, s2) in enumerate(zip(str1.split(';'), str2.split(';'))):
        if s1 != s2:
            return i
    return 8


def shear_results(taxonomy_table: str, taxonomy_map: str, output: str) -> None:
    with open(taxonomy_table, 'r') as inf:
        csv_inf = csv.reader(inf, delimiter="\t")
        columns = next(csv_inf)
        columns = dict(zip(columns[1:], range(len(columns))))
        indptr = [0]
        indices = np.array([], dtype=int)
        data = np.array([], dtype=int)
        names = []
        for ix, row in enumerate(csv_inf):
            if ix % 1000 == 0:
                print(ix)
            names.append(row[0])
            np_row = np.array(row[1:], dtype=int)
            temp_indx = np_row > 0
            data = np.concatenate((data, np_row[temp_indx]))
            indices = np.concatenate((indices, np.where(temp_indx)[0]))
            indptr.append(indices.shape[0])

    csr = csr_matrix((data, indices, indptr), dtype=int).T

    with open(taxonomy_map) as inf:
        csv_inf = csv.reader(inf, delimiter='\t')
        name2taxonomy = dict(csv_inf)

    cols_tax = [name2taxonomy[name] for name in names]
    rows_tax = [name2taxonomy[_.replace(".", "_", 1)] for _ in sorted(columns, key=columns.get)]

    dat = np.zeros((len(rows_tax), 9), dtype=int)
    for i, row_name in enumerate(rows_tax):
        row = csr.getrow(i)
        for j, indx in enumerate(row.indices):
            dat[i, index_lca(rows_tax[i], cols_tax[indx])] += row.data[j]

    print(str(dat[:, 0].sum()))

    df = pd.DataFrame(dat, index=rows_tax)
    df['sum'] = dat.sum(axis=1)
    df.drop(0, axis=1, inplace=True)
    df.to_csv(output, header=False, sep='\t')

    uniqueness_rate_per_level = np.zeros(8, dtype=float)
    for i in range(0, 8):
        # Take the sum of those columns
        num_hits =  df.iloc[:, i].sum()
        # Total number of possible hits
        total_hits = df['sum'].sum()
        # Uniqueness Rate
        uniqueness_rate_per_level[i] = num_hits/total_hits
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    list(zip(levels, uniqueness_rate_per_level))
    print(uniqueness_rate_per_level.sum())


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    shear_results(args.taxonomy_table, args.taxonomy_map, args.output)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
