#! /usr/bin/env python
# Shears fasta into substrings of length w step d
# usage
# shear_db.py db.fasta 100 50
import argparse
import datetime
import os
import re

from helix.utils import download_txt_url

from ete3 import Tree
import dendropy

URL_ARCHAEA = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar122.tree"
URL_BACTERIA = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120.tree"


def make_arg_parser():
    parser = argparse.ArgumentParser(
        description='This is the commandline interface for shear_db',
    )
    parser.add_argument(
        '-a', '--path_archaea', help='Set the file path of the input archaea newick tree', required=False)
    parser.add_argument(
        '-b', '--path_bacteria', help='Set the file path of the input bacteria newick tree', required=False)
    parser.add_argument(
        '-o', '--output', help='Set the directory path of the output combined sequences (default: cwd)',
        default=os.path.join(os.getcwd(), "combo.tree"))
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s ')
    return parser


def gtdb_format_names(str_tree: str) -> str:
    str_tree = re.sub('d__', 'k__', str_tree)
    str_tree = re.sub('RS_', '', str_tree)
    str_tree = re.sub('GB_', '', str_tree)
    return str_tree


def merge_trees(str_arc_tree: str, str_bac_tree: str):
    str_arc_tree = gtdb_format_names(str_arc_tree)

    dendro_tree_arc = dendropy.Tree.get_from_string(
        str_arc_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)

    max_length = 0
    for edge in dendro_tree_arc.postorder_edge_iter():
        if edge.length:
            if edge.length > max_length:
                max_length = edge.length

    str_bac_tree = gtdb_format_names(str_bac_tree)

    dendro_tree_bac = dendropy.Tree.get_from_string(
        str_bac_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)

    for edge in dendro_tree_bac.postorder_edge_iter():
        if edge.length:
            if edge.length > max_length:
                max_length = edge.length

    ete_tree_arc = Tree(dendro_tree_arc.as_string(schema="newick", suppress_rooting=True), format=1, quoted_node_names=True)
    ete_tree_bac = Tree(dendro_tree_bac.as_string(schema="newick", suppress_rooting=True), format=1, quoted_node_names=True)

    ete_tree = Tree(name='root')

    ete_tree.add_child(ete_tree_arc.get_tree_root(), dist=max_length)
    ete_tree.add_child(ete_tree_bac.get_tree_root(), dist=max_length)

    return ete_tree


def main():
    start_time = datetime.datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.abspath(os.path.join(args.output)))
    os.makedirs(outdir, exist_ok=True)

    if not args.path_archaea:
        str_arc_tree = ''.join(_ for _ in download_txt_url(URL_ARCHAEA))
    else:
        with open(args.path_archaea) as inf:
            str_arc_tree = inf.readlines()[0]

    if not args.path_archaea:
        str_bac_tree = ''.join(_ for _ in download_txt_url(URL_BACTERIA))
    else:
        with open(args.path_bacteria) as inf:
            str_bac_tree = inf.readlines()[0]

    ete_tree = merge_trees(str_arc_tree, str_bac_tree)

    ete_tree.write(outfile=args.output, format=3, quoted_node_names=True)

    print("Execution time: %s" % (datetime.datetime.now() - start_time))


if __name__ == '__main__':
    main()
