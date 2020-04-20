import pytest
import os

from helix.filter_fastq import build_read_filter_set
from helix.filter_fastq import filter_reads
from helix.filter_fastq import read_fastq


@pytest.fixture()
def blast():
    return os.path.join("fixtures", "blast_test.b6")


@pytest.fixture()
def fastq():
    return os.path.join("fixtures", "SRR172902.fastq")


# content of test_tmpdir.py
def test_build_read_filter_set(blast):
    test_set = build_read_filter_set(blast)
    assert len(test_set) == 1000
    assert "SRR172902_2416" in test_set


def test_filter_reads(blast, fastq, tmpdir):
    test_set = build_read_filter_set(blast)
    filter_reads(fastq, tmpdir, test_set)
    records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.fq"))))
    assert records == 978
    filtered_records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.filtered.fq"))))
    assert filtered_records == 22
