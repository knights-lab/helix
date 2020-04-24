import pytest
import os

from helix.filter_fastq import build_read_filter_set
from helix.filter_fastq import filter_reads
from helix.filter_fastq import read_fastq


@pytest.fixture()
def blast_combined():
    return os.path.join("fixtures", "blast_test.b6")


@pytest.fixture()
def fastq_combined():
    return os.path.join("fixtures", "SRR172902.fastq")


@pytest.fixture()
def blast_no_change():
    return os.path.join("fixtures", "scratch", "alignment.burst.best.b6")

@pytest.fixture()
def fastq_no_change():
    return os.path.join("fixtures", "scratch", "Study.ID.1.T.0.S1.R1.001.fastq")

# content of test_tmpdir.py
def test_build_read_filter_set_combined(blast_combined):
    test_set = build_read_filter_set(blast_combined)
    assert len(test_set) == 1000
    assert "SRR172902_2416" in test_set


def test_filter_reads_combined(blast_combined, fastq_combined, tmpdir):
    test_set = build_read_filter_set(blast_combined)
    filter_reads(fastq_combined, tmpdir, test_set)
    records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.fq"))))
    assert records == 978
    filtered_records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.filtered.fq"))))
    assert filtered_records == 22


def test_filter_reads_no_change(blast_no_change, fastq_no_change, tmpdir):
    test_set = build_read_filter_set(blast_no_change)
    filter_reads(fastq_no_change, tmpdir, test_set, change_flag=True)
    records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.fq"))))
    assert records == 978
    filtered_records = sum(1 for record in read_fastq(open(os.path.join(tmpdir, "SRR172902.filtered.fq"))))
    assert filtered_records == 22
