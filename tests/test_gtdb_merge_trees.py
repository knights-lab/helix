import pytest
import os
from pathlib import Path

from ete3 import Tree
import dendropy

from helix.gtdb_merge_trees import read_tree, gtdb_format_names


@pytest.fixture()
def tree_archaea() -> Tree:
    return Path("fixtures") / Path("ar122.tree")


@pytest.fixture()
def tree_bacteria() -> Tree:
    return Path("fixtures") / Path("bac120.tree")


def test_combine_trees(tmpdir: Path, tree_archaea: Path, tree_bacteria: Path):
    with open(tree_archaea) as inf:
        str_arc_tree = inf.readlines()[0]

    str_arc_tree = gtdb_format_names(str_arc_tree)

    dendro_tree_arc = dendropy.Tree.get_from_string(
        str_arc_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)

    max_length = 0
    for edge in dendro_tree_arc.postorder_edge_iter():
        if edge.length:
            if edge.length > max_length:
                max_length = edge.length

    with open(tree_bacteria) as inf:
        str_bac_tree = inf.readlines()[0]

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

    ete_tree.write(format=3, outfile=Path(tmpdir) / Path("combo.tree"))

    print('hello, world')

    # terminal_node.add_child(tree_bac)

    # tree_both.get_common_ancestor('GCA_003663295.1', 'GCF_014472415.1')

    # print(dendro_tree_arc)
    # print(tree_archaea)
    # print(tree_bacteria)
