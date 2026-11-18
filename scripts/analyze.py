"""Считает Tm всеми методами для ридов одного организма и сохраняет сводку в CSV.

Пример:
  python scripts/analyze.py --fastq data/ecoli.fastq.gz --name ecoli --n-reads 100000
Для каждого режима (целые риды и окна по W нт) пишется строка на метод:
organism, mode, method, n, gc_mean, len_mean, len_min, len_max, tm_mean, tm_sd.
"""
import argparse
import csv
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.analysis import collect_reads, gc_values, sample_windows, tm_values  # noqa: E402


def summarize(name, mode, seqs, na, ct):
    gcs = gc_values(seqs)
    lens = [len(s) for s in seqs]
    rows = []
    for mname, vals in tm_values(seqs, na=na, ct=ct).items():
        rows.append({
            "organism": name, "mode": mode, "method": mname, "n": len(seqs),
            "gc_mean": round(statistics.fmean(gcs), 3),
            "len_mean": round(statistics.fmean(lens), 2),
            "len_min": min(lens), "len_max": max(lens),
            "tm_mean": round(statistics.fmean(vals), 3),
            "tm_sd": round(statistics.pstdev(vals), 3),
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--name", required=True, help="короткое имя организма, например saureus")
    ap.add_argument("--n-reads", type=int, default=100000)
    ap.add_argument("--min-q", type=float, default=20.0, help="минимальное среднее качество Phred рида")
    ap.add_argument("--windows", type=int, nargs="*", default=[20, 40])
    ap.add_argument("--na", type=float, default=50.0, help="[Na+], мМ")
    ap.add_argument("--ct", type=float, default=50.0, help="суммарная концентрация цепей, нМ")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()

    reads, scanned = collect_reads(a.fastq, a.n_reads, a.min_q)
    print(f"Просмотрено ридов: {scanned}; прошло фильтр: {len(reads)}")
    if not reads:
        sys.exit("Нет ридов, прошедших фильтр (проверьте --min-q и формат файла)")

    lens = [len(r) for r in reads]
    print(f"Длина ридов: средняя {statistics.fmean(lens):.1f}, мин {min(lens)}, макс {max(lens)}")
    rows = summarize(a.name, "full", reads, a.na, a.ct)
    rnd = random.Random(a.seed)  # один генератор на все окна — числа совпадают с прежними результатами
    for w in a.windows:
        rows += summarize(a.name, f"w{w}", sample_windows(reads, w, rnd=rnd), a.na, a.ct)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"summary_{a.name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"Сводка записана: {path}")


if __name__ == "__main__":
    main()
