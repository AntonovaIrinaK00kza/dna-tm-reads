"""Методы расчёта температуры плавления ДНК (Tm, °C).

Все функции принимают последовательность из символов A, C, G, T.
ВНИМАНИЕ: формулы реализованы по публикациям "по памяти" — перед включением в доклад
сверьте их с Biopython (tests/test_tm.py) и первоисточниками.
"""
import math

R = 1.987  # газовая постоянная, кал/(моль*К)

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def gc_percent(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


# ---------- Эмпирические GC-формулы ----------

def tm_wallace(seq):
    """Правило Уоллеса: 2*(A+T) + 4*(G+C). Рассчитано на олигонуклеотиды < ~14 нт."""
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))


def tm_marmur_doty(seq):
    """Marmur–Doty (1962): 69.3 + 0.41*%GC (длинная ДНК, ~0.2 М Na+)."""
    return 69.3 + 0.41 * gc_percent(seq)


def tm_chester_marshak(seq):
    """Формула (1) из ВКР: 69.3 + 0.41*%GC - 650/n."""
    return 69.3 + 0.41 * gc_percent(seq) - 650.0 / len(seq)


def tm_wetmur(seq, na_mM=50.0):
    """Wetmur (1991), ДНК:ДНК: 81.5 + 16.6*log10(Na/(1+0.7*Na)) + 0.41*%GC - 500/n."""
    na = na_mM / 1000.0
    return 81.5 + 16.6 * math.log10(na / (1.0 + 0.7 * na)) + 0.41 * gc_percent(seq) - 500.0 / len(seq)


# ---------- Формула (2) из ВКР (в духе Khandelwal–Bhyravabhotla 2010) ----------

STRENGTH = {
    "GC": 13, "CC": 11, "GG": 11, "CG": 10, "AC": 10, "GT": 10, "TC": 8, "GA": 8,
    "CT": 8, "AG": 8, "TG": 7, "CA": 7, "AT": 7, "TT": 5, "AA": 5, "TA": 4,
}


def tm_khandelwal_like(seq):
    """Tm = 7.31*E + 16.47*ln(n) - 41.32, где E = (сумма значений динуклеотидов)/n.
    Деление на n (а не на n-1) сохранено как в ВКР, чтобы воспроизвести её примеры."""
    n = len(seq)
    total = sum(STRENGTH[seq[i - 1:i + 1]] for i in range(1, n))
    return 7.31 * (total / n) + 16.47 * math.log(n) - 41.32


# ---------- Nearest-neighbor, SantaLucia 1998 (унифицированные параметры) ----------

_NN = {  # (dH ккал/моль, dS кал/(моль*К)), 5'->3'
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}


def tm_nn_santalucia(seq, na_mM=50.0, ct_nM=50.0):
    """Модель ближайших соседей, SantaLucia 1998, с солевой поправкой на Na+.

    ct_nM — суммарная концентрация цепей (по умолчанию 50 нМ, как в Biopython: dnac1=dnac2=25).
    """
    n = len(seq)
    dh = ds = 0.0
    for i in range(n - 1):
        h, s = _NN[seq[i:i + 2]]
        dh += h
        ds += s
    for end in (seq[0], seq[-1]):  # инициация по каждому концу
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1
    self_comp = seq == revcomp(seq)
    if self_comp:
        ds += -1.4
    ds += 0.368 * (n - 1) * math.log(na_mM / 1000.0)
    x = 1.0 if self_comp else 4.0
    ct = ct_nM * 1e-9
    return dh * 1000.0 / (ds + R * math.log(ct / x)) - 273.15


METHODS = {
    "Wallace": tm_wallace,
    "Marmur-Doty": tm_marmur_doty,
    "Chester-Marshak": tm_chester_marshak,
    "Wetmur": tm_wetmur,
    "Khandelwal-like": tm_khandelwal_like,
    "NN SantaLucia 1998": tm_nn_santalucia,
}
