import pytest
from pathlib import Path

from helix.gtdb_merge_trees import merge_trees


@pytest.fixture()
def tree_archaea() -> Path:
    return Path("fixtures") / Path("ar122.tree")


@pytest.fixture()
def tree_bacteria() -> Path:
    return Path("fixtures") / Path("bac120.tree")


def test_combine_trees(tmpdir: Path, tree_archaea: Path, tree_bacteria: Path):
    with open(tree_archaea) as inf:
        str_arc_tree = inf.readlines()[0]

    with open(tree_bacteria) as inf:
        str_bac_tree = inf.readlines()[0]

    ete_tree = merge_trees(str_arc_tree, str_bac_tree)

    ete_tree.write(tmpdir / Path("combo.tree"), format=3)
