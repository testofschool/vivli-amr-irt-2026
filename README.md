# IRT-Based Antimicrobial Resistance Ranking from Sparse Surveillance Data

**2026 Vivli AMR Surveillance Data Challenge** | AMR ID: 00013266

Team: Jung Min Kang, Do Weon Go, Jeong Kyu Choi

## Overview

Item Response Theory (IRT) and the Linear Logistic Test Model (LLTM) applied to the ATLAS Antibiotics dataset (1,011,168 isolates, 83 countries, 2004-2024) via the Vivli AMR Register.

## Key Results

- IRT taxon-level resistance rankings correlate with WHO BPPL 2024 (rho = 0.748, p = 0.020)
- LLTM recovers IRT antibiotic potency ranking (rho = 0.936) using 8 CARD molecular features
- IRT maintains higher ranking stability than averaging under 20-40% data dropout
- Novel drugs (cefiderocol, meropenem-vaborbactam) correctly identified as high-potency

## Repository Structure

```
src/fit_rasch_irt.py              Core IRT + LLTM estimation pipeline
mapping/who_bppl_species_mapping.csv   WHO BPPL 2024 species mapping (n=9)
mapping/card_feature_dictionary.csv    CARD feature extraction schema
figures/                          Final publication figures
outputs/                          Aggregate results (no patient-level data)
```

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
