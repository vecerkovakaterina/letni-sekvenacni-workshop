# Nanopore Summer School Workshop — Environment Setup

This replaces the old, messy conda env / pip requirements. Everything needed
to run `summer_school_workshop.ipynb` is now in `environment.yml`.

## One-time setup

```bash
conda env create -f environment.yml
```

## Every time you want to run the notebook

```bash
conda activate nanopore-workshop
jupyter lab
```

(`conda deactivate` when you're done.)

## Not covered by conda — install separately

- **Dorado** (ONT basecaller used in `print_basecalled.py` and cell 11):
  not on bioconda, download the binary from Oxford Nanopore and put it on
  your `PATH`. Also grab the basecalling model referenced in the notebook
  (`dna_r10.4.1_e8.2_400bps_sup@...`) ahead of time — it's a few GB.
- **IGV Desktop** (used in the "namapovaná čtení" demo, cell ~26): install
  separately and open it manually during that part of the workshop.
- Reference genome FASTA (`Homo_sapiens.GRCh38.dna.toplevel.fa`) and the
  example pod5/fastq/bam files referenced throughout — these are pulled
  from lab paths (`/var/lib/minknow/...`) that need updating to wherever
  this year's run data lands. Search the notebook for `TODO upravit cestu`
  to find every path that needs changing before the session.

## If the env breaks mid-workshop

```bash
conda env remove -n nanopore-workshop
conda env create -f environment.yml
```
