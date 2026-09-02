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
from itertools import combinations

import numpy as np
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from scipy.stats import f_oneway, ttest_ind


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


def plot_chromosome_profile_heatmap(coverages_df_normalized):
    """
    Heatmapa vzorek x chromozom, kde barva = jak moc se normalizované
    pokrytí daného vzorku na daném chromozomu liší od průměru přes
    všechny vzorky (z-skóre).

    Proč tohle a ne korelace mezi okny na jednom chromozomu: celochromo-
    zomální zisk/ztráta posune stejnoměrně ÚROVEŇ pokrytí všech oken na
    tom chromozomu, ale nemění "tvar" (pořadí oken podle pokrytí) uvnitř
    chromozomu - a korelace mezi okny umí zachytit jen tvar, ne úroveň.
    Proto korelace spočítaná zvlášť pro jeden chromozom u celochromozo-
    málních změn nefunguje - je potřeba porovnat přímo ÚROVNĚ pokrytí
    mezi vzorky, což dělá tahle funkce. Navíc díky z-skóre (odečtení
    průměru přes vzorky u každého chromozomu) se automaticky odstraní
    jakýkoli sdílený "genomický" vzor (mapovatelnost, GC obsah), který by
    jinak byl stejný pro všechny vzorky - zůstanou jen skutečné rozdíly
    mezi vzorky/skupinami.

    coverages_df_normalized: výstup z normalize_and_average() (první prvek
        vráceného tuple), tedy pokrytí normalizované podle hloubky
        sekvenování, ale JEŠTĚ ne zprůměrované přes skupiny.

    Vrací: DataFrame (vzorek x chromozom) se z-skóre.
    """
    value_cols = [c for c in coverages_df_normalized.columns if c.endswith("_coverage")]
    sample_names = [c.replace("_coverage", "") for c in value_cols]

    matrix = coverages_df_normalized.set_index("#rname")[value_cols].T
    matrix.index = sample_names
    matrix.columns.name = "chromosome"

    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0)
    z = z.fillna(0.0)

    fig = px.imshow(
        z,
        color_continuous_scale="RdBu_r",
        zmin=-z.abs().max().max(), zmax=z.abs().max().max(),
        aspect="auto",
        title="Profil pokrytí po chromozomech (z-skóre vůči průměru přes všechny vzorky)",
    )
    fig.update_layout(xaxis_title="Chromozom", yaxis_title="Vzorek")
    fig.show()

    return z


MAIN_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def compute_windowed_coverage(bams, reference_genome, read_counts, window=5_000_000, multiply_by=1_000, chromosomes=MAIN_CHROMS):
    """
    Spočítá normalizované pokrytí v oknech (výchozí velikost 5 Mb) podél
    každého chromozomu - to umožní podívat se na pokrytí "zblízka", ne jen
    za celý chromozom.

    Ve výchozím nastavení se počítá jen pro 24 hlavních chromozomů
    (1-22, X, Y) - drobné nezařazené kontigy (typu "KI270539.1") jsou
    jen pár set bází dlouhé, takže do nich prakticky nikdy nic nenamapuje
    a jen by zanášely tabulku samými nulami.

    Velikost okna (`window`) je potřeba přizpůsobit hloubce sekvenování:
    s malým počtem čtení (typické pro krátký workshopový běh) je jemné
    dělení (např. 50 kb) zbytečné - většina oken vyjde nulová bez ohledu
    na to, co je uvnitř. Čím méně čtení máte, tím větší okno zvolte
    (řádově Mb); s hlubším sekvenováním můžete jít na jemnější rozlišení.

    Vrací: depths_df_normalized (chrom, window_start, window_end, sloupec pokrytí pro každý vzorek)
    """
    _run(f"samtools faidx {reference_genome}")
    _run(f"cut -f1,2 {reference_genome}.fai > chrom.sizes")

    sizes = pd.read_csv("chrom.sizes", sep="\t", header=None, names=["chrom", "size"], dtype={"chrom": str})
    if chromosomes is not None:
        sizes = sizes[sizes["chrom"].isin(chromosomes)]
        if sizes.empty:
            raise ValueError(
                "Žádný z požadovaných chromozomů nebyl v referenčním genomu nalezen - "
                "zkontroluj, jak jsou chromozomy pojmenované (např. '1' vs 'chr1')."
            )
    sizes.to_csv("chrom.sizes", sep="\t", header=False, index=False)

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


def plot_sample_similarity(depths_df_normalized):
    """
    Slepá detektivka: spočítá, jak moc jsou si jednotlivé vzorky navzájem
    podobné podle profilu normalizovaného pokrytí napříč všemi chromozomy
    (bez ohledu na to, jestli známe jejich skutečný původ), a vykreslí to
    jako maticový graf (heatmapu).

    Vzorky, které patří ke stejné buněčné linii, by měly mít podobný
    "otisk" chromozomálních aberací, a tedy vysokou korelaci - naopak
    vzorky z různých linií by měly korelovat méně.

    depths_df_normalized: výstup z compute_windowed_coverage()

    Vrací: correlation_df - čtvercová tabulka korelací mezi vzorky
    """
    value_cols = [c for c in depths_df_normalized.columns if c.startswith("overlapping_features_")]
    sample_names = [c.replace("overlapping_features_", "") for c in value_cols]

    correlation_df = depths_df_normalized[value_cols].corr(method="spearman")
    correlation_df.index = sample_names
    correlation_df.columns = sample_names

    fig = px.imshow(
        correlation_df,
        text_auto=".2f",
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        title="Podobnost profilu pokrytí mezi vzorky (Spearmanova korelace)",
    )
    fig.update_layout(xaxis_title="Vzorek", yaxis_title="Vzorek")
    fig.show()

    return correlation_df


def plot_chromosome_coverage(chrom_num, df):
    """Vykreslí normalizované pokrytí jednoho chromozomu pro všechny vzorky."""
    value_cols = [c for c in df.columns if c.startswith("overlapping_features_")]
    df_subset = df.loc[df["chrom"] == str(chrom_num), ]
    fig = px.line(df_subset, x="window_start", y=value_cols, title=f"Pokrytí chromozomu {chrom_num}", line_shape="hv")
    fig.update_layout(xaxis_title="Souřadnice", yaxis_title="Pokrytí")
    fig.show(renderer="notebook")


def run_anova(coverages_df, samples, alpha=0.05):
    """
    Pro každý chromozom porovná pokrytí mezi skupinami vzorků (ANOVA) a vrátí
    chromozomy, kde se skupiny statisticky významně liší (p < alpha).

    DŮLEŽITÉ OMEZENÍ ANOVA: F-test řekne jen "NĚKTERÁ skupina se od ostatních
    liší", ale neřekne KTERÁ a JAKÝM SMĚREM. Pokud chceš zjistit třeba
    "chromozom 17 je zvýšený jen ve skupině U2OS, ne u zdravé kontroly ani
    U251", ANOVA na tohle sama o sobě nestačí - použij run_pairwise_posthoc()
    níže, který porovná skupiny po dvojicích. Proto tahle funkce ke každému
    chromozomu rovnou přidává i průměr pokrytí v každé skupině (sloupce
    mean_<skupina>) - alespoň pro rychlou orientaci, i bez post-hoc testu.

    Navíc kontroluje podezřelý případ: pokud má některá skupina nulový
    rozptyl (tj. všechny repliky mají pro daný chromozom úplně stejnou
    hodnotu pokrytí), vyjde F-statistika "inf" a p-hodnota 0 - ale
    nejde o reálný biologický efekt, jen o dělení nulou. Nejčastější
    příčina je, že dvě "repliky" ve skutečnosti ukazují na stejný FASTQ
    soubor (např. při testování před tím, než jsou k dispozici opravdová
    data). Tyhle řádky se označí ve sloupci `zero_variance_warning` a
    nepočítají se do `significant_df`.

    Vrací: (results_df, significant_df, suspect_df)
      results_df   - všechny chromozomy, včetně mean_<skupina> sloupců
      significant_df - chromozomy s p < alpha A bez podezření na nulový rozptyl
      suspect_df   - chromozomy s p < alpha, ale kde je nulový rozptyl uvnitř
                     některé skupiny (podívej se na tyhle ručně, nejspíš jde o chybu ve vstupních datech)
    """
    groups = {}
    for name, info in samples.items():
        groups.setdefault(info["group"], []).append(f"{name}_coverage")

    results = []
    for _, row in coverages_df.iterrows():
        chrom = row["#rname"]
        group_values = {g: row[cols].values.astype(float) for g, cols in groups.items()}
        zero_variance = any(len(vals) > 1 and np.ptp(vals) == 0 for vals in group_values.values())
        try:
            stat, pvalue = f_oneway(*group_values.values())
        except ZeroDivisionError:
            stat, pvalue = np.inf, 0.0
        row_result = {
            "chromosome": chrom,
            "statistic": stat,
            "pvalue": pvalue,
            "zero_variance_warning": zero_variance,
        }
        for g, vals in group_values.items():
            row_result[f"mean_{g}"] = vals.mean()
        results.append(row_result)

    results_df = pd.DataFrame(results).sort_values("pvalue")
    below_alpha = results_df["pvalue"] < alpha
    significant_df = results_df[below_alpha & ~results_df["zero_variance_warning"]]
    suspect_df = results_df[below_alpha & results_df["zero_variance_warning"]]

    if len(suspect_df) > 0:
        print(
            f"⚠️  {len(suspect_df)} chromozomů má p < {alpha}, ale s podezřením na "
            "nulový rozptyl uvnitř skupiny (viz suspect_df) - zkontroluj, jestli "
            "repliky ve `samples` opravdu ukazují na různé FASTQ soubory."
        )

    return results_df, significant_df, suspect_df


def run_pairwise_posthoc(coverages_df, samples, chromosomes=None, alpha=0.05):
    """
    Post-hoc test navazující na run_anova(): pro zadané chromozomy porovná
    KAŽDOU DVOJICI skupin zvlášť (Welchův t-test, nepředpokládá stejný
    rozptyl mezi skupinami) - to řekne, KTERÁ konkrétní skupina se od které
    liší a jakým směrem. Např. pro "chromozom 17 zvýšený jen v jedné
    skupině" bys tu měl/a vidět něco jako:

        chromosome=17, group_a=cancer_U2OS, group_b=healthy_control,
        direction="cancer_U2OS > healthy_control", pvalue_bonf < 0.05
        chromosome=17, group_a=cancer_U2OS, group_b=cancer_U251,
        direction="cancer_U2OS > cancer_U251", pvalue_bonf < 0.05
        chromosome=17, group_a=healthy_control, group_b=cancer_U251,
        direction="...", pvalue_bonf >= 0.05  (tahle dvojice se neliší)

    Proč Bonferroniho korekce: pro každý chromozom se dělá víc testů
    najednou (jedna dvojice skupin = jeden test), a čím víc testů, tím
    větší šance na "významný" výsledek jen náhodou. Bonferroniho korekce
    (pvalue_bonf = min(1, pvalue * počet_dvojic_skupin)) tohle kompenzuje -
    je to konzervativní (spíš podhodnotí významnost), ale jednoduchá a
    bezpečná volba.

    ⚠️ Se 2 replikami na skupinu má i tenhle test malou statistickou sílu -
    ber výsledky spíš jako orientační vodítko, kam se dál dívat (např. do
    IGV), než jako definitivní důkaz.

    coverages_df: výstup z compute_genome_coverage() (NEnormalizovaný -
        stejně jako u run_anova - pro post-hoc porovnání skupin mezi sebou
        na tom nezáleží, normalizace posouvá všechny skupiny stejně).
    chromosomes: které chromozomy testovat (default: všechny v coverages_df -
        v praxi má smysl sem dát jen chromozomy z significant_df, aby se
        zbytečně nenafukoval počet testů).

    Vrací: DataFrame se sloupci chromosome, group_a, group_b, mean_a, mean_b,
    direction, pvalue, pvalue_bonf, significant_bonf - seřazený podle pvalue.
    """
    groups = {}
    for name, info in samples.items():
        groups.setdefault(info["group"], []).append(f"{name}_coverage")

    if chromosomes is None:
        subset_df = coverages_df
    else:
        wanted = {str(c) for c in chromosomes}
        subset_df = coverages_df[coverages_df["#rname"].astype(str).isin(wanted)]

    group_pairs = list(combinations(sorted(groups), 2))
    n_comparisons = max(1, len(group_pairs))

    results = []
    for _, row in subset_df.iterrows():
        chrom = row["#rname"]
        for group_a, group_b in group_pairs:
            values_a = row[groups[group_a]].values.astype(float)
            values_b = row[groups[group_b]].values.astype(float)
            mean_a, mean_b = values_a.mean(), values_b.mean()

            if np.ptp(values_a) == 0 and np.ptp(values_b) == 0:
                # obě skupiny mají nulový rozptyl - t-test by dělil nulou;
                # pokud jsou navíc stejné, není mezi nimi žádný rozdíl k testování
                stat, pvalue = (0.0, 1.0) if mean_a == mean_b else (np.inf, 0.0)
            else:
                stat, pvalue = ttest_ind(values_a, values_b, equal_var=False)

            pvalue_bonf = min(1.0, pvalue * n_comparisons)
            if mean_a > mean_b:
                direction = f"{group_a} > {group_b}"
            elif mean_b > mean_a:
                direction = f"{group_b} > {group_a}"
            else:
                direction = f"{group_a} = {group_b}"

            results.append({
                "chromosome": chrom,
                "group_a": group_a,
                "group_b": group_b,
                "mean_a": mean_a,
                "mean_b": mean_b,
                "direction": direction,
                "pvalue": pvalue,
                "pvalue_bonf": pvalue_bonf,
                "significant_bonf": pvalue_bonf < alpha,
            })

    return pd.DataFrame(results).sort_values("pvalue")
