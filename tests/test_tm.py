import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools import tm_methods as tm  # noqa: E402

EX1 = "A" * 26
EX2 = "AACGCGCGTGAATTCTGGCCGCAA"


def test_vkr_example_1_and_2_formula_2():
    # значения из ВКР: 47.48 и 73.76 (в тексте ВКР для примера 2 ошибочно напечатано 74.6)
    assert tm.tm_khandelwal_like(EX1) == pytest.approx(47.48, abs=0.02)
    assert tm.tm_khandelwal_like(EX2) == pytest.approx(73.76, abs=0.02)


def test_gc_formulas():
    seq = "GC" * 10  # 100 % GC, n = 20
    assert tm.tm_marmur_doty(seq) == pytest.approx(69.3 + 41.0)
    assert tm.tm_chester_marshak(seq) == pytest.approx(69.3 + 41.0 - 650 / 20)
    assert tm.tm_wallace("ATGC") == 2 * 2 + 4 * 2


def test_nn_is_monotonic_in_gc():
    low = "AT" * 10 + "AA"
    high = "GC" * 10 + "GG"
    assert tm.tm_nn_santalucia(low) < tm.tm_nn_santalucia(high)


def test_nn_matches_biopython_if_installed():
    mt = pytest.importorskip("Bio.SeqUtils.MeltingTemp")
    seqs = ["AACGTGAATTCTGGCA", EX2, "GCGATCGCATGCAAGTCCGA", "ATATATATTTAAAATATTAT"]
    for s in seqs:
        ref = mt.Tm_NN(s, nn_table=mt.DNA_NN4, Na=50, dnac1=25, dnac2=25, saltcorr=5)
        mine = tm.tm_nn_santalucia(s, na_mM=50, ct_nM=50)
        # допуск 1 °C: возможны различия в обработке концов; большее расхождение — повод искать ошибку
        assert abs(ref - mine) < 1.0, (s, ref, mine)
