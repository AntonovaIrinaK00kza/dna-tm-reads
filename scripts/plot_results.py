"""Строит графики «Tm от GC» по results/summary_*.csv: по одной панели на режим (full, w20, w40).

Примеры:
  python scripts/plot_results.py
  python scripts/plot_results.py --exclude Wallace --modes w20 w40 --out tm_vs_gc_clean.png
Метод можно убрать, если его значения нефизичны и сжимают шкалу остальных (например, Wallace на длинных ридах).
Линии соединяют точки для наглядности, это не регрессия; планки — SD между ридами, а не ошибка среднего.
"""

import argparse
import csv
import glob
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default="results")
    ap.add_argument(
        "--exclude", nargs="*", default=[], help="методы, которые не рисовать"
    )
    ap.add_argument(
        "--modes",
        nargs="*",
        default=None,
        help="какие режимы рисовать, например full w40",
    )
    ap.add_argument("--out", default="tm_vs_gc.png")
    a = ap.parse_args()

    data = defaultdict(
        lambda: defaultdict(list)
    )  # mode -> method -> [(gc, tm, sd, organism)]
    for p in sorted(glob.glob(str(Path(a.results) / "summary_*.csv"))):
        with open(p, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["method"] in a.exclude:
                    continue
                data[r["mode"]][r["method"]].append(
                    (
                        float(r["gc_mean"]),
                        float(r["tm_mean"]),
                        float(r["tm_sd"]),
                        r["organism"],
                    )
                )
    if not data:
        sys.exit(
            "Не найдено results/summary_*.csv — сначала запустите scripts/analyze.py"
        )

    modes = sorted(data, key=lambda m: (m != "full", m))
    if a.modes:
        modes = [m for m in modes if m in a.modes]
    fig, axes = plt.subplots(
        1, len(modes), figsize=(5 * len(modes), 4.2), squeeze=False
    )
    # разные маркеры и штрихи: сборник печатается, и рисунок должен читаться в чёрно-белом виде
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    styles = ["-", "--", "-.", ":", "-", "--", "-."]
    order = sorted({m for mode in data for m in data[mode]})
    for ax, mode in zip(axes[0], modes):
        for method, pts in data[mode].items():
            pts = sorted(pts)
            k = order.index(method)
            ax.errorbar(
                [x[0] for x in pts],
                [x[1] for x in pts],
                yerr=[x[2] for x in pts],
                marker=markers[k % len(markers)],
                linestyle=styles[k % len(styles)],
                capsize=2,
                label=method,
                color=f"C{k}",
            )
        ax.set_title("Целые риды" if mode == "full" else f"Окна {mode[1:]} нт")
        ax.set_xlabel("Среднее GC, %")
        ax.set_ylabel("Tm, °C")
        ax.grid(alpha=0.3)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=min(len(labels), 5),
        fontsize=8,
        frameon=False,
    )
    fig.tight_layout(
        rect=(0, 0.08, 1, 1)
    )  # легенда под панелями, не перекрывает данные
    out = Path(a.results) / a.out
    fig.savefig(out, dpi=150)
    print(f"Сохранено: {out}")


if __name__ == "__main__":
    main()
