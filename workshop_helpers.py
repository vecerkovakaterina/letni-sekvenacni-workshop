"""
Pomocné funkce pro workshop o Nanopore sekvenování.

Tady je "pod kapotou" veškerý bioinformatický kód (spojování FASTQ souborů,
zarovnání na referenční genom pomocí minimap2, počítání pokrytí pomocí
samtools/bedtools, normalizace, statistika). V notebooku studenti volají
jen tyhle funkce a dívají se na výsledky - nemusí číst desítky řádků
shellových příkazů.
"""

import glob
import subprocess

import numpy as np
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from scipy.stats import f_oneway


def _run(command):
    """Spustí shellový příkaz a vyhodí chybu se srozumitelnou hláškou, pokud selže."""
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Příkaz selhal: {command}\n{result.stderr}")


def merge_and_align(samples, reference_genome, threads=12):
    """
    Pro každý vzorek: spojí všechny FASTQ soubory do jednoho, zarovná je na
    referenční genom (minimap2), seřadí a vyfiltruje pouze primárně
    namapovaná čtení (samtools).

    samples: dict {nazev_vzorku: {"fastq_glob": "cesta/k/souborum_*.fastq.gz", "group": "skupina"}}
    reference_genome: cesta k referenčnímu genomu (FASTA, může být .gz)

    Vrací: dict {nazev_vzorku: cesta_k_vyfiltrovanemu_bam_souboru}
    """
    bams = {}
    for name, info in samples.items():
        merged_fastq = f"{name}.fastq.gz"
        sorted_bam = f"{name}.sorted.bam"
        filtered_bam = f"{name}.primary_mapped.sorted.bam"

        matched_files = glob.glob(info["fastq_glob"])
        if not matched_files:
            raise FileNotFoundError(f"Nenalezeny žádné soubory pro vzorek '{name}' podle vzoru {info['fastq_glob']}")

        _run(f"cat {' '.join(matched_files)} > {merged_fastq}")
        _run(f"minimap2 -ax map-ont -t {threads} \"{reference_genome}\" {merged_fastq} | samtools sort -o {sorted_bam} && samtools index {sorted_bam}")
        _run(f"samtools view -F 2308 {sorted_bam} -o {filtered_bam} && samtools index {filtered_bam}")

        bams[name] = filtered_bam
        print(f"✓ {name}: zarovnáno a vyfiltrováno")

    return bams


def _count_reads(fastq_path):
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))


def compute_genome_coverage(bams):
    """
    Pro každý vzorek spočítá pokrytí jednotlivých chromozomů (samtools coverage)
    a počet namapovaných čtení (potřebný pro normalizaci).

    bams: dict {nazev_vzorku: cesta_k_bam_souboru}, výstup z merge_and_align()

    Vrací: (coverages_df, read_counts)
      coverages_df - tabulka s pokrytím každého chromozomu pro každý vzorek
      read_counts  - dict {nazev_vzorku: pocet_namapovanych_cteni}
    """
    per_sample_dfs = []
    read_counts = {}

    for name, bam in bams.items():
        coverage_tsv = f"{name}.coverage.tsv"
        _run(f"samtools coverage {bam} | head -n 25 > {coverage_tsv}")

        df = pd.read_csv(coverage_tsv, sep="\t", usecols=["#rname", "startpos", "endpos", "numreads"])
        df[f"{name}_coverage"] = df["numreads"] / df["endpos"] # normalizovane na delku chromozomu
        df = df.drop(columns="numreads")
        per_sample_dfs.append(df)

        fastq_path = f"{name}.primary_mapped.fastq"
        _run(f"samtools fastq {bam} > {fastq_path}")
        read_counts[name] = _count_reads(fastq_path)

    coverages_df = per_sample_dfs[0]
    for df in per_sample_dfs[1:]:
        coverages_df = coverages_df.merge(df)

    return coverages_df, read_counts


def normalize_and_average(coverages_df, read_counts, samples, multiply_by=1_000_000_000):
    """
    Normalizuje pokrytí podle počtu namapovaných čtení v každém vzorku (aby šly
    vzorky s různou hloubkou sekvenování spravedlivě porovnat) a spočítá
    průměr pokrytí přes replikáty v rámci každé skupiny (zdravý/nádor 1/nádor 2...).

    Vrací: (coverages_df_normalized, coverages_mean_df)
    """
    coverages_df_normalized = coverages_df.copy(deep=True)
    for name, count in read_counts.items():
        col = f"{name}_coverage"
        coverages_df_normalized[col] = (coverages_df[col] / count) * multiply_by

    groups = {}
    for name, info in samples.items():
        groups.setdefault(info["group"], []).append(f"{name}_coverage")

    coverages_mean_df = coverages_df_normalized[["#rname"]].copy()
    for group, cols in groups.items():
        coverages_mean_df[f"{group}_mean"] = coverages_df_normalized[cols].mean(axis=1)

    return coverages_df_normalized, coverages_mean_df


def plot_genome_overview(coverages_mean_df):
    """Vykreslí průměrné (normalizované) pokrytí všech chromozomů pro každou skupinu vzorků."""
    group_cols = [c for c in coverages_mean_df.columns if c != "#rname"]
    fig = px.scatter(coverages_mean_df, x="#rname", y=group_cols, title="Pokrytí chromozomů zprůměrované přes replikáty")
    fig.update_traces(marker_size=10)
    fig.update_layout(xaxis_title="Chromozom #", yaxis_title="Pokrytí")
    fig.show()


def compute_windowed_coverage(bams, reference_genome, read_counts, window=50_000, multiply_by=1_000):
    """
    Spočítá normalizované pokrytí v oknech (např. po 50 kb) podél každého
    chromozomu - to umožní podívat se na pokrytí "zblízka", ne jen za celý
    chromozom.

    Vrací: depths_df_normalized (chrom, window_start, window_end, sloupec pokrytí pro každý vzorek)
    """
    _run(f"samtools faidx {reference_genome}")
    _run(f"cut -f1,2 {reference_genome}.fai > chrom.sizes")
    _run(f"bedtools makewindows -g chrom.sizes -w {window} > {window}bps_windows.bed")

    per_sample_dfs = []
    for name, bam in bams.items():
        windowed_tsv = f"{name}.coverage_{window}bps_windows.tsv"
        _run(f"bedtools coverage -a {window}bps_windows.bed -b {bam} > {windowed_tsv}")
        df = pd.read_csv(
            windowed_tsv, header=None, sep="\t", low_memory=False,
            usecols=[0, 1, 2, 3],
            names=["chrom", "window_start", "window_end", f"overlapping_features_{name}"],
        )
        per_sample_dfs.append(df)

    merge_cols = ["chrom", "window_start", "window_end"]
    depths_df = per_sample_dfs[0]
    for df in per_sample_dfs[1:]:
        depths_df = depths_df.merge(df, on=merge_cols)

    depths_df_normalized = depths_df.copy(deep=True)
    for name, count in read_counts.items():
        col = f"overlapping_features_{name}"
        depths_df_normalized[col] = (depths_df[col] / count) * multiply_by

    return depths_df_normalized


def plot_chromosome_coverage(chrom_num, df):
    """Vykreslí normalizované pokrytí jednoho chromozomu pro všechny vzorky."""
    value_cols = [c for c in df.columns if c.startswith("overlapping_features_")]
    df_subset = df.loc[df["chrom"] == str(chrom_num), ]
    fig = px.line(df_subset, x="window_start", y=value_cols, title=f"Pokrytí chromozomu {chrom_num}", line_shape="hv")
    fig.update_layout(xaxis_title="Souřadnice", yaxis_title="Pokrytí")
    fig.show()


def run_anova(coverages_df, samples, alpha=0.05):
    """
    Pro každý chromozom porovná pokrytí mezi skupinami vzorků (ANOVA) a vrátí
    chromozomy, kde se skupiny statisticky významně liší (p < alpha).

    Vrací: (results_df, significant_df)
    """
    groups = {}
    for name, info in samples.items():
        groups.setdefault(info["group"], []).append(f"{name}_coverage")

    results = []
    for _, row in coverages_df.iterrows():
        chrom = row["#rname"]
        group_values = [row[cols].values for cols in groups.values()]
        try:
            stat, pvalue = f_oneway(*group_values)
        except ZeroDivisionError:
            stat, pvalue = np.inf, 0.00001
        results.append({"chromosome": chrom, "statistic": stat, "pvalue": pvalue})

    results_df = pd.DataFrame(results).sort_values("pvalue")
    significant_df = results_df[results_df["pvalue"] < alpha]
    return results_df, significant_df
