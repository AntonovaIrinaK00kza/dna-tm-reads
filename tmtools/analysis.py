"""Общая логика анализа ридов: используется и командной строкой (scripts/analyze.py), и Streamlit-приложением (app.py)."""
import random
import statistics

from .fastq import is_clean, mean_phred, read_fastq
from .tm_methods import METHODS, gc_percent, tm_nn_santalucia, tm_wetmur


def collect_reads(path, n_reads, min_q=20.0):
    """Первые n_reads ридов, прошедших фильтр (только ACGT, средний Phred >= min_q). Возвращает (риды, просмотрено)."""
    reads, scanned = [], 0
    for _, seq, qual in read_fastq(path):
        scanned += 1
        if is_clean(seq) and mean_phred(qual) >= min_q:
            reads.append(seq)
            if len(reads) >= n_reads:
                break
    return reads, scanned


def sample_windows(reads, w, seed=1, rnd=None):
    """По одному случайному окну длины w из каждого рида (риды короче w пропускаются).

    rnd — общий генератор random.Random: analyze.py передаёт один и тот же объект для всех размеров окон,
    чтобы воспроизвести прежние числа; без него генератор создаётся из seed."""
    if rnd is None:
        rnd = random.Random(seed)
    out = []
    for r in reads:
        if len(r) >= w:
            p = rnd.randrange(0, len(r) - w + 1)
            out.append(r[p:p + w])
    return out


def gc_values(seqs):
    return [gc_percent(s) for s in seqs]


def tm_values(seqs, na=50.0, ct=50.0, methods=None):
    """Словарь {метод: [Tm для каждой последовательности]} в условиях na (мМ) и ct (нМ)."""
    result = {}
    for name, fn in METHODS.items():
        if methods is not None and name not in methods:
            continue
        if fn is tm_nn_santalucia:
            result[name] = [fn(s, na_mM=na, ct_nM=ct) for s in seqs]
        elif fn is tm_wetmur:
            result[name] = [fn(s, na_mM=na) for s in seqs]
        else:
            result[name] = [fn(s) for s in seqs]
    return result


def describe(values):
    """Среднее, SD (генеральное, как в analyze.py), минимум, максимум."""
    return statistics.fmean(values), statistics.pstdev(values), min(values), max(values)
