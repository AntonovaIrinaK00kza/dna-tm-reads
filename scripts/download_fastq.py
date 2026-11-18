"""Скачивает только НАЧАЛО файла .fastq.gz по ссылке и сохраняет его как обычный .fastq.
Работает на Windows, Linux и macOS без gunzip/head/curl.

Пример:
  python scripts/download_fastq.py "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR270/003/SRR27019103/SRR27019103_1.fastq.gz" --reads 200000 --out data/saureus.fastq
"""
import argparse
import gzip
import sys
import urllib.request


def normalize(url):
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if "://" not in url:
        return "https://" + url
    return url


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url", help="ссылка на *_1.fastq.gz (из колонки FASTQ FTP в ENA)")
    ap.add_argument("--reads", type=int, default=200000, help="сколько ридов взять (по 4 строки на рид)")
    ap.add_argument("--out", required=True, help="куда сохранить, например data/saureus.fastq")
    a = ap.parse_args()

    url = normalize(a.url.strip().strip('"'))
    max_lines = a.reads * 4
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    written = 0
    with urllib.request.urlopen(req, timeout=60) as resp, \
            gzip.open(resp, "rt", encoding="utf-8") as src, \
            open(a.out, "w", encoding="utf-8", newline="\n") as dst:
        for line in src:
            if written >= max_lines:
                break
            dst.write(line)
            written += 1
    print(f"Записано ридов: {written // 4} -> {a.out}")
    if written % 4:
        sys.exit("Внимание: число строк не кратно 4, файл мог оборваться")


if __name__ == "__main__":
    main()
