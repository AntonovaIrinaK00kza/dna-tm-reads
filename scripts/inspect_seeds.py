"""Диагностика k-мерного поиска: какие k-меры порождают найденные пары.

Если несколько k-меров дают основную часть пар (повторы, рРНК, плазмиды) или это низкосложные k-меры
(например, из одних A/T), число пар нельзя трактовать как «честные перекрытия ридов в геноме».

Пример: python scripts/inspect_seeds.py --fastq data/saureus.fastq --n 5000 --k 20 --step 10
"""

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tmtools.analysis import collect_reads  # noqa: E402
from tmtools.tm_methods import revcomp  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--min-q", type=float, default=20.0)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument(
        "--dump",
        default="largest_group_reads.fasta",
        help="куда записать несколько ридов самой большой группы (для BLAST)",
    )
    a = ap.parse_args()

    reads, _ = collect_reads(a.fastq, a.n, a.min_q)
    index = defaultdict(set)
    for i, r in enumerate(reads):
        for p in range(0, len(r) - a.k + 1, a.step):
            index[r[p : p + a.k]].add(i)

    pair_kmers = defaultdict(
        set
    )  # пара ридов -> канонические k-меры, которые её породили
    for i, r in enumerate(reads):
        rc = revcomp(r)
        for p in range(len(rc) - a.k + 1):
            km = rc[p : p + a.k]
            for j in index.get(km, ()):
                if j != i:
                    canon = min(km, revcomp(km))
                    pair_kmers[(min(i, j), max(i, j))].add(canon)

    total = len(pair_kmers)
    print(f"Ридов: {len(reads)}; пар найдено: {total}")
    if not total:
        return
    support = Counter()
    for kms in pair_kmers.values():
        for km in kms:
            support[km] += 1
    low = {km for km in support if len(set(km)) <= 2}
    only_low = sum(1 for kms in pair_kmers.values() if kms <= low)
    print(
        f"Пар, объяснимых ТОЛЬКО низкосложными k-мерами (<= 2 разных нуклеотида): {only_low} ({100 * only_low / total:.1f} %)"
    )
    top = [km for km, _ in support.most_common(a.top)]
    covered = sum(1 for kms in pair_kmers.values() if kms & set(top))
    print(
        f"Пар, в которых есть хотя бы один из {a.top} самых частых k-меров: {covered} ({100 * covered / total:.1f} %)"
    )
    # группировка пар по «участкам»: связные компоненты графа пар (union-find)
    parent = list(range(len(reads)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in pair_kmers:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    involved = {x for pair in pair_kmers for x in pair}
    comp_reads, comp_pairs = Counter(), Counter()
    for x in involved:
        comp_reads[find(x)] += 1
    for i, _ in pair_kmers:
        comp_pairs[find(i)] += 1
    print(
        f"\nРидов, участвующих хотя бы в одной паре: {len(involved)} из {len(reads)} ({100 * len(involved) / len(reads):.1f} %)"
    )
    print(f"Групп связанных ридов (вероятных участков генома): {len(comp_reads)}")
    print("Самые большие группы (ридов / пар):")
    for root, n_reads in comp_reads.most_common(5):
        print(
            f"  {n_reads} ридов, {comp_pairs[root]} пар ({100 * comp_pairs[root] / total:.1f} % всех пар)"
        )
    # сравнение самой большой группы с остальными ридами и выгрузка ридов для BLAST
    big_root = comp_reads.most_common(1)[0][0]
    big = [i for i in involved if find(i) == big_root]
    rest = [i for i in range(len(reads)) if find(i) != big_root or i not in involved]
    gc = (
        lambda seqs: 100
        * sum(x.count("G") + x.count("C") for x in seqs)
        / sum(len(x) for x in seqs)
    )
    print(
        f"\nСредний GC ридов самой большой группы: {gc([reads[i] for i in big]):.1f} %; "
        f"остальных ридов: {gc([reads[i] for i in rest]):.1f} %"
    )
    with open(a.dump, "w", encoding="utf-8") as f:
        for i in sorted(big)[:5]:
            f.write(f">read_{i}\n{reads[i]}\n")
    print(
        f"5 ридов самой большой группы записаны в {a.dump}: вставьте их в NCBI BLAST (blastn)"
    )
    print(f"\nТоп-{a.top} k-меров по числу пар:")
    for km, c in support.most_common(a.top):
        print(f"  {km}  пар: {c}  (разных нуклеотидов: {len(set(km))})")


if __name__ == "__main__":
    main()
