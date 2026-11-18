"""Сравнивает NN-реализацию (SantaLucia 1998) с Biopython на случайных последовательностях
и печатает максимальное и среднее расхождение — число для доклада.

Пример: python scripts/compare_biopython.py --n 1000 --lengths 16 20 30 40
"""
import argparse
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.tm_methods import tm_nn_santalucia  # noqa: E402

try:
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    sys.exit("Biopython не установлен: pip install biopython")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000, help="число случайных последовательностей на длину")
    ap.add_argument("--lengths", type=int, nargs="*", default=[16, 20, 30, 40])
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    rnd = random.Random(a.seed)
    for length in a.lengths:
        diffs = []
        worst = (0.0, "")
        for _ in range(a.n):
            s = "".join(rnd.choices("ACGT", k=length))
            ref = mt.Tm_NN(s, nn_table=mt.DNA_NN4, Na=a.na, dnac1=25, dnac2=25, saltcorr=5)
            mine = tm_nn_santalucia(s, na_mM=a.na, ct_nM=50.0)
            d = abs(ref - mine)
            diffs.append(d)
            if d > worst[0]:
                worst = (d, s)
        print(f"длина {length}: n={a.n}, макс. расхождение {max(diffs):.3f} °C, "
              f"среднее {statistics.fmean(diffs):.3f} °C (худший случай: {worst[1]})")


if __name__ == "__main__":
    main()
