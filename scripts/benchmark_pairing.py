"""Сравнение алгоритмов поиска комплементарных ридов: перебор, хэширование, сортировка, k-мерные «семена».

Замеряется время (и, с --memory, пиковая память по tracemalloc).
Примеры:
  python scripts/benchmark_pairing.py --max-brute 20000 --max-n 200000            # синтетические риды
  python scripts/benchmark_pairing.py --fastq data/saureus.fastq --max-n 20000     # реальные риды
  python scripts/benchmark_pairing.py --max-brute 4000 --max-n 200000 --memory
"""
import argparse
import csv
import random
import sys
import time
import tracemalloc
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.analysis import collect_reads  # noqa: E402
from tmtools.pairing import brute_force_pairs, hash_pairs, kmer_seed_pairs, sort_pairs  # noqa: E402
from tmtools.tm_methods import revcomp  # noqa: E402

SIZES = (500, 1000, 2000, 4000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)


def make_reads(n, length, planted, seed=1):
    rnd = random.Random(seed)
    reads = ["".join(rnd.choices("ACGT", k=length)) for _ in range(n)]
    for _ in range(int(n * planted)):  # подсаживаем известные комплементарные пары
        i = rnd.randrange(n)
        j = rnd.randrange(n)
        reads[j] = revcomp(reads[i])
    return reads


def timed(fn, reads):
    t = time.perf_counter()
    res = fn(reads)
    return time.perf_counter() - t, res


def peak_mib(fn, reads):
    """Пиковая дополнительная память (МиБ), выделенная вызовом fn; риды должны быть созданы до вызова."""
    tracemalloc.start()
    base, _ = tracemalloc.get_traced_memory()
    tracemalloc.reset_peak()
    fn(reads)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return (peak - base) / 2 ** 20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", help="взять риды из файла вместо синтетических")
    ap.add_argument("--min-q", type=float, default=20.0)
    ap.add_argument("--length", type=int, default=100, help="длина синтетических ридов")
    ap.add_argument("--max-brute", type=int, default=4000)
    ap.add_argument("--max-seed", type=int, default=50000, help="максимум N для k-мерного поиска")
    ap.add_argument("--max-n", type=int, default=200000)
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--memory", action="store_true", help="измерить пиковую память (медленнее)")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()

    if a.fastq:
        pool, _ = collect_reads(a.fastq, a.max_n, a.min_q)
        sizes = [n for n in SIZES if n <= len(pool)] or [len(pool)]
        source = "real"
    else:
        pool, sizes, source = None, [n for n in SIZES if n <= a.max_n], "synthetic"

    rows = []
    for n in sizes:
        reads = pool[:n] if pool is not None else make_reads(n, a.length, planted=0.01)
        t_hash, p_hash = timed(hash_pairs, reads)
        t_sort, p_sort = timed(sort_pairs, reads)
        assert p_hash == p_sort, "hash и sort дали разные результаты"
        t_brute = None
        if n <= a.max_brute:
            t_brute, p_brute = timed(brute_force_pairs, reads)
            assert p_brute == p_hash, "перебор и hash дали разные результаты"
        t_seed = n_seed = None
        if n <= a.max_seed:
            t_seed, p_seed = timed(lambda r: kmer_seed_pairs(r, k=a.k, step=a.step), reads)
            n_seed = len(p_seed)
        row = {"source": source, "n": n, "brute_s": None if t_brute is None else round(t_brute, 4),
               "hash_s": round(t_hash, 4), "sort_s": round(t_sort, 4),
               "seed_s": None if t_seed is None else round(t_seed, 4),
               "pairs_exact": len(p_hash), "pairs_seed": n_seed}
        if a.memory:
            row["mem_hash_MiB"] = round(peak_mib(hash_pairs, reads), 2)
            row["mem_sort_MiB"] = round(peak_mib(sort_pairs, reads), 2)
            row["mem_brute_MiB"] = round(peak_mib(brute_force_pairs, reads), 2) if t_brute is not None else None
            row["mem_seed_MiB"] = round(peak_mib(lambda r: kmer_seed_pairs(r, k=a.k, step=a.step), reads), 2) \
                if t_seed is not None else None
        rows.append(row)
        print(row)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "benchmark_pairing.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)

    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    series = [("brute_s", "Перебор O(N²)", "o-"), ("hash_s", "Хэширование O(N)", "s--"),
              ("sort_s", "Сортировка O(N log N)", "^-."), ("seed_s", f"k-мерные семена (k={a.k}, шаг {a.step})", "d:")]
    for key, label, style in series:
        pts = [(r["n"], r[key]) for r in rows if r[key] is not None]
        if pts:
            ax.loglog([p[0] for p in pts], [p[1] for p in pts], style, label=label)
    ax.set_xlabel("Число ридов N")
    ax.set_ylabel("Время, с")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "benchmark_pairing.png", dpi=150)
    print("Сохранено:", out / "benchmark_pairing.png")


if __name__ == "__main__":
    main()
