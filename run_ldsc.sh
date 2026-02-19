#!/usr/bin/env bash
set -euo pipefail

# Activate conda environment
conda activate ldsc

# Munge ADHD sumstats
./munge_sumstats.py \
  --sumstats adhd_excEAF.qc \
  --out adhd_demontis_22_munge \
  --merge-alleles w_hm3.snplist \
  --chunksize 500000

# Munge telomere sumstats
./munge_sumstats.py \
  --sumstats telomere_excEAF.qc \
  --N 472174 \
  --out telomere_munge \
  --merge-alleles w_hm3.snplist \
  --chunksize 500000

# Run LDSC genetic correlation
./ldsc.py \
  --rg telomere_munge.sumstats.gz,adhd_demontis_22_munge.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out Telomere21_ADHD22
