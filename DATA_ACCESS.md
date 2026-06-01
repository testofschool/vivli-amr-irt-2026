# Data Access

## ATLAS Antibiotics (Pfizer via Vivli)

Obtained through Vivli AMR Register (https://amr.vivli.org), Data Request ID 00013266, approved May 2026.

**This dataset cannot be redistributed.** To reproduce:

1. Create account at https://amr.vivli.org
2. Request ATLAS_Antibiotics dataset
3. After approval, download the CSV
4. Place at data/atlas_vivli_2004_2024.csv
5. Run: python src/fit_rasch_irt.py

Our May 2026 download: 1,011,168 isolates, 83 countries, 2004-2024.
Public Vivli metadata may list a different count from an earlier snapshot.

## CARD v4.0.0

Freely available: https://card.mcmaster.ca/download

```bash
curl -sL https://card.mcmaster.ca/latest/data -o card_data.tar.bz2
tar xjf card_data.tar.bz2
```

Our feature extraction yielded 6,402 antibiotic-linked resistance gene entries from aro_categories_index.tsv.

## Global AMR R&D Hub

Dashboard export (May 2026): https://dashboard.globalamrhub.org
