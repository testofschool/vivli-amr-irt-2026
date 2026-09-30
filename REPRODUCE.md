# Reproduction Guide

## What this repository contains

- `src/fit_rasch_irt.py` — CLI pipeline: fits 1PL IRT + LLTM to an ATLAS matrix
- `mapping/` — WHO BPPL 2024 species mapping and CARD feature dictionary
- `figures/` — final publication figures
- `outputs/` — **derived aggregate results** (already included; no patient-level data)

## What you must supply yourself

### 1. ATLAS extract (NOT redistributable)

Raw ATLAS data cannot be shared under the Vivli data use terms. To obtain it:

1. Register at https://amr.vivli.org
2. Submit a data request for `ATLAS_Antibiotics`
3. After approval, download the CSV extract

The analysis in this repo used the **Vivli-approved May 2026 release**
(1,011,168 isolates, 83 countries, 2004–2024). Older public Vivli metadata
snapshots may list a different isolate count.

### 2. CARD ontology

Freely available:

```bash
curl -sL https://card.mcmaster.ca/latest/data -o card.tar.bz2
tar xjf card.tar.bz2
# key file: aro_categories_index.tsv
```

CARD's official release contains far more than the 6,402 entries we use; 6,402
is the **parsed antibiotic-linked subset** produced by our feature extraction,
not CARD's total size.

## Running the pipeline

```bash
pip install -r requirements.txt

python src/fit_rasch_irt.py \
    --atlas /path/to/atlas_extract.csv \
    --card  /path/to/aro_categories_index.tsv \
    --out   outputs/
```

The script exits gracefully with an explanatory message if either input is missing.

## Outputs

| File | Description |
|------|-------------|
| `outputs/superbug_leaderboard.csv` | Species ranked by latent resistance theta |
| `outputs/antibiotic_potency.csv`   | Antibiotics ranked by potency b (IRT + LLTM) |

> **Note (2026-09-29):** the files already committed in `outputs/` (`ATLAS_*.csv`) and `figures/`
> do not match what this script writes, and the code and command that produced them are not in
> this repository (different file names/columns; committed θ and b are not mean-centred; the
> sparsity, temporal and geographic outputs have no committed generator). See README
> "Reproducibility status of committed outputs".

## Method notes

- **Likelihood:** each observed species–antibiotic cell contributes
  `r_ij ~ Binomial(n_ij, p_ij)`, with `n_ij` capped at 500 during weighting so
  high-volume cells do not dominate estimation.
- **Identifiability:** theta and b are mean-centered after optimization
  (82 effective degrees of freedom for a 41+42 model).
- **Validation:** IRT theta vs WHO BPPL 2024 priority tiers, Spearman rho = 0.748
  (p = 0.020, n = 9); dropping the imperfect *E. faecalis* proxy gives rho = 0.718
  (p = 0.045, n = 8). See `mapping/who_bppl_species_mapping.csv`.
