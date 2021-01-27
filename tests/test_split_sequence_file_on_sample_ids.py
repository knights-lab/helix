import pytest
import os

from helix.split_sequence_file_on_sample_ids import split_sequence_file


@pytest.fixture
def fasta() -> str:
    return os.path.join("fixtures", "combined_seqs.200k.fna")


def test_split_sequence_file(fasta, tmpdir):
    buffer = 10_000
    split_sequence_file(fasta, tmpdir, buffer)
