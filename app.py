"""Streamlit-приложение: распределения GC и Tm по ридам FASTQ и сравнение алгоритмов поиска комплементарных ридов.

Запуск (из папки проекта):  streamlit run app.py
"""
import tempfile
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from tmtools.analysis import collect_reads, describe, gc_values, sample_windows, tm_values  # noqa: E402
from tmtools.pairing import brute_force_pairs, hash_pairs, kmer_seed_pairs, sort_pairs  # noqa: E402
from tmtools.tm_methods import METHODS  # noqa: E402

st.set_page_config(page_title="Tm ДНК по ридам", layout="wide")
st.title("Температура плавления ДНК по ридам секвенирования")
st.caption("Распределения GC и Tm по шести методам; сравнение алгоритмов поиска комплементарных ридов.")


@st.cache_data(show_spinner=False)
def load_reads(path, n_reads, min_q):
    return collect_reads(path, n_reads, min_q)


@st.cache_data(show_spinner=False)
def compute_tm(path, n_reads, min_q, mode_w, na, ct, methods):
    """mode_w = 0 — целые риды, иначе длина окна. Кэш по параметрам: повторный запуск с теми же значениями мгновенный."""
    reads, _ = collect_reads(path, n_reads, min_q)
    seqs = reads if mode_w == 0 else sample_windows(reads, mode_w)
    return gc_values(seqs), tm_values(seqs, na=na, ct=ct, methods=list(methods)), [len(s) for s in seqs]


def timed_best(fn, data, repeats=3):
    """Лучшее время из нескольких запусков (для быстрых алгоритмов): снижает шум от прогрева и сборщика мусора."""
    best, res = None, None
    for _ in range(repeats):
        t0 = time.perf_counter()
        res = fn(data)
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
        if dt > 1.0:  # медленные запуски не повторяем
            break
    return best, res


def resolve_path(uploaded, path_text):
    """Загруженный файл сохраняем во временный (один раз на файл); иначе берём путь с диска."""
    if uploaded is not None:
        key = (uploaded.name, uploaded.size)
        if st.session_state.get("upload_key") != key:
            suffix = ".fastq.gz" if uploaded.name.endswith(".gz") else ".fastq"
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
            tmp.write(uploaded.getbuffer())
            tmp.close()
            st.session_state["upload_key"] = key
            st.session_state["upload_path"] = tmp.name
        return st.session_state["upload_path"]
    if path_text.strip():
        return path_text.strip().strip('"')
    return None


# ---------------- боковая панель ----------------
st.sidebar.header("Данные")
uploaded = st.sidebar.file_uploader("FASTQ или FASTQ.GZ", type=["fastq", "fq", "gz"])
path_text = st.sidebar.text_input("…или путь к файлу на диске", value="")
n_reads = st.sidebar.slider("Сколько ридов взять", 1000, 200000, 20000, step=1000)
min_q = st.sidebar.slider("Минимальный средний Phred", 0, 40, 20)

st.sidebar.header("Условия расчёта Tm")
na = st.sidebar.number_input("[Na⁺], мМ", min_value=1.0, max_value=1000.0, value=50.0)
ct = st.sidebar.number_input("Концентрация цепей, нМ", min_value=1.0, max_value=100000.0, value=50.0)

path = resolve_path(uploaded, path_text)
if path is None:
    st.info("Загрузите файл FASTQ слева или укажите путь к нему.")
    st.stop()
if not Path(path).is_file():
    st.error(f"Файл не найден: {path}")
    st.stop()

with st.spinner("Читаю FASTQ…"):
    reads, scanned = load_reads(path, n_reads, min_q)
if not reads:
    st.error("Не нашлось ридов, прошедших фильтр: проверьте формат файла и порог качества.")
    st.stop()

mean_len = sum(len(r) for r in reads) / len(reads)
c1, c2, c3 = st.columns(3)
c1.metric("Ридов после фильтра", f"{len(reads):,}".replace(",", " "))
c2.metric("Просмотрено ридов", f"{scanned:,}".replace(",", " "))
c3.metric("Средняя длина, нт", f"{mean_len:.1f}")

tab_tm, tab_pairs, tab_about = st.tabs(["Распределения Tm и GC", "Поиск комплементарных ридов", "О проекте"])

# ---------------- вкладка 1: распределения ----------------
with tab_tm:
    left, right = st.columns([1, 2])
    with left:
        mode = st.radio("Что считать", ["Целые риды", "Окно фиксированной длины"])
        window = st.slider("Длина окна, нт", 14, 100, 40) if mode != "Целые риды" else 0
        names = list(METHODS)
        chosen = st.multiselect("Методы", names, default=[n for n in names if n != "Wallace"])
        st.caption("Правило Уоллеса по умолчанию скрыто: на длинных последовательностях оно даёт нефизичные значения.")
    if not chosen:
        st.warning("Выберите хотя бы один метод.")
    else:
        with st.spinner("Считаю Tm…"):
            gcs, tms, lens = compute_tm(path, n_reads, min_q, window, na, ct, tuple(chosen))
        with right:
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.hist(gcs, bins=50, color="C0", alpha=0.85)
            ax.set_xlabel("GC, %")
            ax.set_ylabel("Число последовательностей")
            ax.set_title("Распределение GC")
            fig.tight_layout()
            st.pyplot(fig)
        fig, ax = plt.subplots(figsize=(9, 3.8))
        for i, (name, vals) in enumerate(tms.items()):
            ax.hist(vals, bins=60, alpha=0.5, label=name, color=f"C{i}")
        ax.set_xlabel("Tm, °C")
        ax.set_ylabel("Число последовательностей")
        ax.set_title("Распределения Tm по методам")
        ax.legend(fontsize=8)
        fig.tight_layout()
        st.pyplot(fig)

        rows = []
        for name, vals in tms.items():
            mean, sd, lo, hi = describe(vals)
            rows.append({"Метод": name, "Среднее, °C": round(mean, 2), "SD, °C": round(sd, 2),
                         "Мин, °C": round(lo, 2), "Макс, °C": round(hi, 2)})
        table = pd.DataFrame(rows)
        st.dataframe(table, use_container_width=True)
        st.download_button("Скачать таблицу (CSV)", table.to_csv(index=False).encode("utf-8"),
                           file_name="tm_summary.csv", mime="text/csv")

# ---------------- вкладка 2: алгоритмы ----------------
with tab_pairs:
    st.write("Два рида образуют дуплекс, если один равен обратному комплементу другого. "
             "Сравниваем время и число найденных пар для разных алгоритмов на ридах из загруженного файла.")
    p1, p2, p3 = st.columns(3)
    n_pairs = p1.slider("Число ридов N", 500, len(reads), min(len(reads), 5000), step=500) if len(reads) > 500 else len(reads)
    k = p2.slider("Длина k-мера", 10, 40, 20)
    step = p3.slider("Шаг прореживания", 1, 20, 10)
    include_brute = st.checkbox("Включить полный перебор O(N²) (медленно при N > 5000)", value=n_pairs <= 5000)
    if st.button("Запустить сравнение"):
        subset = reads[:n_pairs]
        algos = [("Хэширование", "O(N)", hash_pairs), ("Сортировка", "O(N log N)", sort_pairs),
                 (f"k-мерные семена (k={k}, шаг {step})", "≈ O(N·L/шаг)", lambda r: kmer_seed_pairs(r, k=k, step=step))]
        if include_brute:
            algos.insert(0, ("Полный перебор", "O(N²)", brute_force_pairs))
        results, found = [], {}
        for name, cx, fn in algos:
            with st.spinner(f"{name}…"):
                dt, res = timed_best(fn, subset)
            found[name] = res
            results.append({"Алгоритм": name, "Сложность": cx, "Время, с": round(dt, 4), "Найдено пар": len(res)})
        st.dataframe(pd.DataFrame(results), use_container_width=True)
        exact = [found[n] for n in ("Полный перебор", "Хэширование", "Сортировка") if n in found]
        if all(e == exact[0] for e in exact):
            st.success("Точные алгоритмы нашли одинаковые наборы пар.")
        else:
            st.error("Точные алгоритмы дали разные наборы пар: это ошибка реализации.")
        fig, ax = plt.subplots(figsize=(7, 3.2))
        ax.barh([r["Алгоритм"] for r in results], [r["Время, с"] for r in results])
        ax.set_xscale("log")
        ax.set_xlabel("Время, с (лог. шкала)")
        fig.tight_layout()
        st.pyplot(fig)
        st.caption("В реальных данных полностью комплементарных ридов мало; частичную комплементарность "
                   "(общий участок) находит k-мерный поиск — он приближённый и тратит больше памяти. "
                   "Время — лучшее из нескольких запусков (медленные не повторяются).")

# ---------------- вкладка 3: о проекте ----------------
with tab_about:
    st.markdown(
        """
**Методы Tm:** Marmur–Doty, Chester–Marshak, Wetmur, формула из предшествующей работы (Khandelwal-like),
nearest-neighbor (SantaLucia 1998) и правило Уоллеса. Условия расчёта задаются слева.

**Ограничения:** экспериментальных Tm нет, поэтому методы сравниваются между собой, а не с реальностью;
формулы для олигонуклеотидов на длинных ридах выходят за область применимости.

**Код и тесты:** пакет `tmtools`, тесты `pytest`, сверка NN-реализации с Biopython (`scripts/compare_biopython.py`).
"""
    )
