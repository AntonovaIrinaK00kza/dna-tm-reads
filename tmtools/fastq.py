"""Потоковое чтение FASTQ.

Важно: запись FASTQ — это ровно 4 строки (заголовок, последовательность, '+', качество).
Последовательность определяется по ПОЗИЦИИ строки, а не по первой букве: строки качества
тоже могут начинаться с A, C, G или T (это допустимые символы Phred).
"""
import gzip


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def read_fastq(path):
    """Генератор (заголовок, последовательность, качество). Память не растёт с размером файла."""
    with open_text(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip("\r\n")
            plus = f.readline()
            qual = f.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Нарушен формат FASTQ около записи: {header[:40]!r}")
            if len(seq) != len(qual):
                raise ValueError(f"Длины последовательности и качества не совпадают: {header[:40]!r}")
            yield header[1:].rstrip("\r\n"), seq.upper(), qual


def mean_phred(qual, offset=33):
    """Среднее качество Phred (кодировка Phred+33, как у современного Illumina)."""
    if not qual:
        return 0.0
    return sum(ord(c) - offset for c in qual) / len(qual)


def is_clean(seq):
    """Только A, C, G, T (без N и прочих символов)."""
    return bool(seq) and set(seq) <= set("ACGT")
