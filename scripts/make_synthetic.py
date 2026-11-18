"""Генерирует тестовый FASTQ с заданным GC-составом (чтобы отладить конвейер без скачивания данных).

Пример: python scripts/make_synthetic.py --gc 0.33 --reads 20000 --out data/synth_low.fastq
"""
import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.tm_methods import revcomp  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gc", type=float, default=0.5, help="доля GC в геноме (0..1)")
    ap.add_argument("--reads", type=int, default=20000)
    ap.add_argument("--length", type=int, default=150)
    ap.add_argument("--genome", type=int, default=200000, help="длина синтетического генома")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rnd = random.Random(a.seed)
    w = [(1 - a.gc) / 2, a.gc / 2, a.gc / 2, (1 - a.gc) / 2]  # A, C, G, T
    genome = "".join(rnd.choices("ACGT", weights=w, k=a.genome))
    with open(a.out, "w", encoding="utf-8") as f:
        for i in range(a.reads):
            p = rnd.randrange(0, a.genome - a.length)
            read = genome[p:p + a.length]
            if rnd.random() < 0.5:
                read = revcomp(read)
            # качество: в основном высокое, иногда низкое; строка может начинаться с A/C/G/T
            qual = "".join(rnd.choice("AAAABBFFFF") for _ in range(a.length)) if i % 10 else "#" * a.length
            f.write(f"@synthetic_{i}\n{read}\n+\n{qual}\n")
    print(f"Записано {a.reads} ридов в {a.out}")


if __name__ == "__main__":
    main()
