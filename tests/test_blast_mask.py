import os

import pytest
import numpy as np

from helix.blast_mask import gen_blast_reads_hits, build_mask_dict, mask_fasta, read_fasta


@pytest.fixture
def gtdb_20_database() -> str:
    return os.path.join("fixtures", "gtdb_20.fna")


@pytest.fixture
def gtdb_20_blast() -> str:
    return os.path.join("fixtures", "gtdb_20.1k.b6")


def test_read_gtdb_20(gtdb_20_blast):
    gen = gen_blast_reads_hits(gtdb_20_blast)
    pass


@pytest.fixture
def gtdb_20_blast_dict(gtdb_20_blast):
    dd = build_mask_dict(gtdb_20_blast)
    return dd


def test_mask_blast(gtdb_20_database, gtdb_20_blast_dict):
    d_sums_before = {}
    with open(gtdb_20_database) as inf:
        gen_fasta = read_fasta(inf)
        for title, data in gen_fasta:
            np_data = np.fromiter(data, dtype='S1')
            d_sums_before[title] = np.sum(np_data == b"N")

    with open(gtdb_20_database) as inf:
        gen_fasta = read_fasta(inf)
        gen_mask_fasta = mask_fasta(gen_fasta, gtdb_20_blast_dict)
        for header, data in gen_mask_fasta:
            np_data = np.fromiter(data, dtype="S1")
            s = np.sum(np_data == b"N")

            num_N_in_masks = 0
            for mask in gtdb_20_blast_dict[header]:
                num_N_in_masks += (mask[1] - mask[0])

            assert d_sums_before[header] + num_N_in_masks == s
