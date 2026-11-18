import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.pairing import brute_force_pairs, hash_pairs, kmer_seed_pairs, sort_pairs  # noqa: E402
from tmtools.tm_methods import revcomp  # noqa: E402


def test_revcomp():
    assert revcomp("AACG") == "CGTT"


def test_three_methods_agree_including_palindromes_and_duplicates():
    rnd = random.Random(3)
    reads = ["".join(rnd.choices("ACGT", k=12)) for _ in range(60)]
    reads += [revcomp(reads[0]), revcomp(reads[1]), reads[0], "ACGT", "ACGT", "GAATTC"]
    b, h, s = brute_force_pairs(reads), hash_pairs(reads), sort_pairs(reads)
    assert b == h == s
    assert (0, 60) in b


def test_position_by_position_complement_is_not_a_pair():
    # 'AAAA' и 'TTTT' — комплементарны и антипараллельны; 'ACGT' с 'TGCA' — не пара (TGCA != revcomp(ACGT)=ACGT)
    assert brute_force_pairs(["AAAA", "TTTT"]) == {(0, 1)}
    assert brute_force_pairs(["ACGT", "TGCA"]) == set()


def test_kmer_seed_finds_partial_overlap():
    a = "ACGTTGCAAGCTTGGATCCA" + "GATTACAGATTACAGATTAC"
    b = revcomp(a[:25]) + "TTTTTTTTTT"
    assert (0, 1) in kmer_seed_pairs([a, b], k=20)
