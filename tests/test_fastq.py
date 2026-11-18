import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.fastq import is_clean, mean_phred, read_fastq  # noqa: E402


def test_quality_line_starting_with_acgt_is_not_a_sequence(tmp_path):
    # строка качества начинается с 'A' — старый фильтр startswith(('A','T','C','G')) принял бы её за рид
    p = tmp_path / "t.fastq"
    p.write_text("@r1\nACGTACGT\n+\nAAABBFFF\n@r2\nTTTTGGGG\n+\nGGGGFFFF\n", encoding="utf-8")
    recs = list(read_fastq(p))
    assert [s for _, s, _ in recs] == ["ACGTACGT", "TTTTGGGG"]


def test_clean_and_phred():
    assert is_clean("ACGT") and not is_clean("ACGN") and not is_clean("")
    assert mean_phred("IIII") == 40
