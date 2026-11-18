import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.analysis import collect_reads, describe, gc_values, sample_windows, tm_values  # noqa: E402
from tmtools.pairing import kmer_seed_pairs  # noqa: E402
from tmtools.tm_methods import METHODS, revcomp  # noqa: E402


def _write(tmp_path, recs):
    p = tmp_path / "t.fastq"
    p.write_text("".join(f"@r{i}\n{s}\n+\n{q}\n" for i, (s, q) in enumerate(recs)), encoding="utf-8")
    return p


def test_collect_reads_filters_n_and_quality(tmp_path):
    p = _write(tmp_path, [("ACGTACGTAC", "I" * 10), ("ACGTNCGTAC", "I" * 10), ("ACGTACGTAC", "#" * 10),
                          ("GGGGCCCCAA", "I" * 10)])
    reads, scanned = collect_reads(p, 10, 20)
    assert reads == ["ACGTACGTAC", "GGGGCCCCAA"] and scanned == 4


def test_collect_reads_limit(tmp_path):
    p = _write(tmp_path, [("ACGTACGTAC", "I" * 10)] * 5)
    reads, scanned = collect_reads(p, 2, 20)
    assert len(reads) == 2 and scanned == 2


def test_sample_windows_length_and_skip_short():
    reads = ["ACGT" * 10, "ACG", "TTTT" * 12]
    w = sample_windows(reads, 20, seed=3)
    assert len(w) == 2 and all(len(x) == 20 for x in w)
    assert sample_windows(reads, 20, seed=3) == w  # воспроизводимость


def test_tm_values_has_all_methods_and_matches_direct_call():
    seqs = ["AACGCGCGTGAATTCTGGCCGCAA", "GATTACAGATTACAGATTAC"]
    res = tm_values(seqs)
    assert set(res) == set(METHODS)
    assert res["Chester-Marshak"][0] == METHODS["Chester-Marshak"](seqs[0])
    only = tm_values(seqs, methods=["Wetmur"])
    assert list(only) == ["Wetmur"]


def test_describe_and_gc_values():
    assert gc_values(["GGCC", "AATT"]) == [100.0, 0.0]
    mean, sd, lo, hi = describe([1.0, 3.0])
    assert (mean, sd, lo, hi) == (2.0, 1.0, 1.0, 3.0)


def test_kmer_step_finds_long_overlap_and_is_subset_of_step1():
    rnd = random.Random(5)
    a = "".join(rnd.choices("ACGT", k=80))
    b = revcomp(a[10:70]) + "".join(rnd.choices("ACGT", k=20))  # общий участок 60 нт
    others = ["".join(rnd.choices("ACGT", k=80)) for _ in range(30)]
    reads = [a, b] + others
    full = kmer_seed_pairs(reads, k=20, step=1)
    sparse = kmer_seed_pairs(reads, k=20, step=10)
    assert (0, 1) in sparse and sparse <= full
