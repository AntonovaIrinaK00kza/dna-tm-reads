"""Поиск полностью комплементарных ридов (антипараллельных: рид A == обратный комплемент рида B).

Три реализации с одинаковым результатом и разной сложностью:
  brute_force_pairs — O(N^2) (как в program1.py, но с правильным обратным комплементом)
  hash_pairs        — O(N) по словарю
  sort_pairs        — O(N log N) по сортировке канонических форм
kmer_seed_pairs — частичная комплементарность (общий k-мер с обратным комплементом), приближённый поиск.
"""
from collections import defaultdict

from .tm_methods import revcomp


def brute_force_pairs(reads):
    rcs = [revcomp(r) for r in reads]
    pairs = set()
    n = len(reads)
    for i in range(n):
        ri = reads[i]
        for j in range(i + 1, n):
            if ri == rcs[j]:
                pairs.add((i, j))
    return pairs


def hash_pairs(reads):
    index = defaultdict(list)
    for i, r in enumerate(reads):
        index[r].append(i)
    pairs = set()
    for i, r in enumerate(reads):
        for j in index.get(revcomp(r), ()):
            if i < j:
                pairs.add((i, j))
    return pairs


def sort_pairs(reads):
    keyed = []
    for i, r in enumerate(reads):
        rc = revcomp(r)
        keyed.append((min(r, rc), i))
    keyed.sort()
    pairs = set()
    start = 0
    n = len(keyed)
    while start < n:
        end = start
        while end < n and keyed[end][0] == keyed[start][0]:
            end += 1
        group = [i for _, i in keyed[start:end]]
        canon = keyed[start][0]
        if revcomp(canon) == canon:  # палиндром: все риды в группе парные
            for a in range(len(group)):
                for b in range(a + 1, len(group)):
                    pairs.add((min(group[a], group[b]), max(group[a], group[b])))
        else:
            side_a = [i for i in group if reads[i] == canon]
            side_b = [i for i in group if reads[i] != canon]
            for a in side_a:
                for b in side_b:
                    pairs.add((min(a, b), max(a, b)))
        start = end
    return pairs


def kmer_seed_pairs(reads, k=20, step=1):
    """Пары ридов с общим комплементарным участком: k-мер одного рида совпадает с k-мером обратного комплемента другого.

    step > 1 — прореживание: в индекс кладутся k-меры только с позиций 0, step, 2*step, ...
    (память и время уменьшаются примерно в step раз); находятся пары с комплементарным участком
    длиной не менее k + step - 1. Это приближённый поиск «семенами», как в выравнивателях.
    """
    index = defaultdict(set)
    for i, r in enumerate(reads):
        for p in range(0, len(r) - k + 1, step):
            index[r[p:p + k]].add(i)
    pairs = set()
    for i, r in enumerate(reads):
        rc = revcomp(r)
        for p in range(len(rc) - k + 1):
            for j in index.get(rc[p:p + k], ()):
                if j != i:
                    pairs.add((min(i, j), max(i, j)))
    return pairs
