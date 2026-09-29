# IRT-Based Antimicrobial Resistance Ranking from Sparse Surveillance Data

**2026 Vivli AMR Surveillance Data Challenge** | AMR ID: 00013266

Team: Jung Min Kang, Do Weon Go, Jeong Kyu Choi

## Overview

Item Response Theory (IRT) and the Linear Logistic Test Model (LLTM) applied to the ATLAS Antibiotics dataset (1,011,168 isolates, 83 countries, 2004-2024) via the Vivli AMR Register.

## Key Results

- IRT taxon-level resistance rankings correlate with WHO BPPL 2024 (rho = 0.748, p = 0.020)
- LLTM recovers IRT antibiotic potency ranking (rho = 0.936) using 8 CARD molecular features
- IRT maintains higher mean ranking stability than averaging at 10-40% additional data dropout (largest gap at 20-40%), but **not** at 50%, where averaging is higher (see Sparsity stress test below)
- Novel drugs (cefiderocol, meropenem-vaborbactam) correctly identified as high-potency

## Repository Structure

```
src/fit_rasch_irt.py              Core IRT + LLTM estimation pipeline
mapping/who_bppl_species_mapping.csv   WHO BPPL 2024 species mapping (n=9)
mapping/card_feature_dictionary.csv    CARD feature extraction schema
figures/                          Final publication figures
outputs/                          Aggregate results (no patient-level data)
```

## Sparsity stress test: exact values and range (corrected 2026-09-29)

From `outputs/ATLAS_sparsity_stress_test.csv` (Spearman rank recovery, mean ± SD over 5 seeds):

| Additional dropout | IRT | Averaging | IRT − Averaging |
|---|---|---|---|
| 10% | 0.982 ± 0.006 | 0.981 ± 0.010 | +0.001 |
| 20% | 0.969 ± 0.010 | 0.952 ± 0.017 | +0.017 |
| 30% | 0.964 ± 0.019 | 0.938 ± 0.019 | +0.026 |
| 40% | 0.940 ± 0.016 | 0.913 ± 0.030 | +0.027 |
| 50% | **0.880** ± 0.014 | **0.916** ± 0.029 | **−0.036** |

IRT degrades more slowly than averaging only from 10% to 40% additional dropout (the 10% gap is
negligible and the SD bands overlap at every level). At 50% dropout the ordering reverses:
IRT 0.8798701298701298 vs averaging 0.9162337662337661 (exact CSV values). The title of
`figures/fig3_sparsity_stress_test.png` ("IRT Degrades More Slowly") and its x-axis stop at 40% and
omit the 50% row; read the figure together with this table. The figure is left unmodified because the
code that drew it is not in this repository.

## Reproducibility status of committed outputs (2026-09-29)

The files in `outputs/` and `figures/` were **not** produced by the committed
`src/fit_rasch_irt.py`; they come from a different, uncommitted version of the analysis code and
**cannot be regenerated from the committed script**. The ATLAS input is not redistributable and is
not in this repository, so the committed script could not be re-run to settle this; the evidence
below comes from the committed files alone:

- `fit_rasch_irt.py` mean-centres θ and b after fitting (`theta -= theta.mean()`, `b -= b.mean()`),
  but in `outputs/ATLAS_superbug_leaderboard.csv` the mean θ over the 41 species is −1.273, and in
  `outputs/ATLAS_drug_potency.csv` the mean b over the 42 antibiotics is +1.243 — neither is centred.
- File names and columns differ: the script writes `superbug_leaderboard.csv`
  (`rank, species, theta, avg_resistance, n_antibiotics`) and `antibiotic_potency.csv`
  (`rank, antibiotic, b_irt, b_lltm, n_species`); the committed files are
  `ATLAS_superbug_leaderboard.csv` (`rank, species, theta, avg_r, n_abx, n_tests`) and
  `ATLAS_drug_potency.csv` (`rank, antibiotic, b, n_spp, R_pct`).
- No committed code produces `ATLAS_sparsity_stress_test.csv`, `ATLAS_temporal_trends.csv`,
  `ATLAS_geographic_trends.csv`, or any of the four figures.

Affected outputs: all five files in `outputs/` and all four files in `figures/`. Because θ and b
are only identified up to a location shift, rankings (and Spearman ρ values computed from them) are
unaffected by the missing centring, but absolute θ/b values in the committed CSVs are on a different
scale from what the committed script would print.

## Data Access

ATLAS data cannot be redistributed. See DATA_ACCESS.md for instructions.
CARD v4.0.0: https://card.mcmaster.ca/download (freely available)

## Requirements

```
pip install -r requirements.txt
```

## Citation

This publication is based on research using data from Pfizer, obtained through https://amr.vivli.org

## Related Work

- Kang (2026). The Scaling Law of Evaluation Failure. arXiv:2605.11205
- Fischer (1973). The Linear Logistic Test Model. Acta Psychologica.

## License

MIT License
