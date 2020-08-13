import os

import pytest

from helix.strip_taxamap import gen_strip_dangling_taxa, strip_dangling_taxa

# @pytest.fixture
# def alltax() -> str:
#     return os.path.join("fixtures", "alltax.txt")
#
#
# def test_strip_dangling_taxa(alltax):
#     gen = gen_strip_dangling_taxa(alltax)
#     for i in gen:
#         print(i)


def test_strip_dangling_taxa_single():
    inp = "k__Archaea;p__Crenarchaeota;c__;o__Acidilobales;f__Caldisphaeraceae;g__Caldisphaera;s__Caldisphaera lagunensis;t__"
    exp = "k__Archaea;p__Crenarchaeota;c__;o__Acidilobales;f__Caldisphaeraceae;g__Caldisphaera;s__Caldisphaera lagunensis"

    results = strip_dangling_taxa(inp)
    assert results == exp


def test_strip_dangling_taxa_double():
    inp = "k__Archaea;p__Crenarchaeota;c__;o__Acidilobales;f__Caldisphaeraceae;g__Caldisphaera;s__;t__"
    exp = "k__Archaea;p__Crenarchaeota;c__;o__Acidilobales;f__Caldisphaeraceae;g__Caldisphaera"

    results = strip_dangling_taxa(inp)
    assert results == exp
