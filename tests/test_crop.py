import os

import pytest

from helix.crop import crop, read_fasta


@pytest.fixture
def gtdb_20_database() -> str:
    return os.path.join("fixtures", "gtdb_20.fna")


def test_crop_200_200(gtdb_20_database, tmpdir):
    outf = os.path.join(tmpdir, "test_crop.fna")

    crop(gtdb_20_database, outf, 200, 200)

    with open(outf) as inf:
        for header, seq in read_fasta(inf):
            assert len(seq) == 200


def test_crop_200_100(gtdb_20_database, tmpdir):
    outf = os.path.join(tmpdir, "test_crop.fna")

    crop(gtdb_20_database, outf, 100, 200)

    with open(outf) as inf:
        for header, seq in read_fasta(inf):
            assert 100 <= len(seq) <= 200
