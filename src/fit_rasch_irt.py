#!/usr/bin/env python3
"""
AMR-IRT: Item Response Theory for Antimicrobial Resistance Ranking
2026 Vivli AMR Surveillance Data Challenge | AMR ID 00013266
Team: Jung Min Kang, Do Weon Go, Jeong Kyu Choi

Fits a 1PL Rasch IRT model and an LLTM (CARD-feature) model to an ATLAS
species x antibiotic resistance matrix.

USAGE:
    python fit_rasch_irt.py --atlas /path/to/atlas_extract.csv \
                            --card  /path/to/card/aro_categories_index.tsv \
                            --out   outputs/

The ATLAS extract is NOT redistributable. Obtain it via the Vivli AMR
Register (https://amr.vivli.org). See REPRODUCE.md for details.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr


def load_atlas_matrix(atlas_path, min_isolates_species=1000, min_tests_cell=30):
    """Build a species x antibiotic resistance matrix from an ATLAS CSV extract."""
    df = pd.read_csv(atlas_path, low_memory=False)
    interp_cols = sorted([c for c in df.columns if c.endswith('_I') and
        not any(c.startswith(g) for g in ['ACC', 'ACT', 'AMPC', 'CMY', 'CTX', 'DHA',
        'FOX', 'GES', 'GIM', 'IMP', 'KPC', 'NDM', 'OXA', 'PER', 'SHV', 'SPM', 'TEM',
        'VEB', 'VIM'])])

    counts = df['Species'].value_counts()
    species = counts[counts >= min_isolates_species].index.tolist()

    records = []
    for col in interp_cols:
        abx = col.replace('_I', '')
        for sp in species:
            vals = df.loc[df['Species'] == sp, col].dropna().astype(str).str.strip().str.lower()
            n = len(vals)
            if n >= min_tests_cell:
                records.append({'species': sp, 'antibiotic': abx,
                                'n': n, 'rate': (vals == 'resistant').sum() / n})
    agg = pd.DataFrame(records)
    pathogens = sorted(agg['species'].unique())
    antibiotics = sorted(agg['antibiotic'].unique())
    pidx = {p: i for i, p in enumerate(pathogens)}
    aidx = {a: i for i, a in enumerate(antibiotics)}

    J, I = len(pathogens), len(antibiotics)
    R = np.full((J, I), np.nan)
    N = np.zeros((J, I))
    for _, r in agg.iterrows():
        R[pidx[r['species']], aidx[r['antibiotic']]] = r['rate']
        N[pidx[r['species']], aidx[r['antibiotic']]] = r['n']
    return pathogens, antibiotics, R, N, ~np.isnan(R)


def fit_irt_1pl(R, mask, N, cap=500, reg=0.01):
    """1PL Rasch: P(R=1) = sigma(theta_j - b_i). Binomial likelihood, n capped at `cap`."""
    J, I = R.shape
    js, is_ = np.where(mask)
    ns = np.minimum(N[js, is_], cap)
    ss = np.round(R[js, is_] * ns)

    def nll(params):
        th, b = params[:J], params[J:]
        p = np.clip(expit(th[js] - b[is_]), 1e-12, 1 - 1e-12)
        return -np.sum(ss * np.log(p) + (ns - ss) * np.log(1 - p)) + reg * np.sum(params**2)

    def grad(params):
        th, b = params[:J], params[J:]
        p = np.clip(expit(th[js] - b[is_]), 1e-12, 1 - 1e-12)
        d = ss - ns * p
        g = np.zeros(J + I)
        np.add.at(g[:J], js, -d)
        np.add.at(g[J:], is_, d)
        return g + 2 * reg * params

    x0 = np.zeros(J + I)
    res = minimize(nll, x0, jac=grad, method='L-BFGS-B',
                   options={'maxiter': 5000, 'ftol': 1e-14})
    theta, b = res.x[:J], res.x[J:]
    # Resolve location non-identifiability: mean-center (82 effective DOF for 41+42)
    theta -= theta.mean()
    b -= b.mean()
    return theta, b, res


def fit_lltm(R, mask, N, X, cap=500, reg=0.01):
    """LLTM: b_i = sum_f w_f * X_if. Joint estimation of theta and feature weights."""
    J, I = R.shape
    F = X.shape[1]
    Xm, Xs = X.mean(0), X.std(0)
    Xs[Xs == 0] = 1
    Xn = (X - Xm) / Xs
    js, is_ = np.where(mask)
    ns = np.minimum(N[js, is_], cap)
    ss = np.round(R[js, is_] * ns)

    def nll(params):
        th, w = params[:J], params[J:J + F]
        b = Xn @ w
        p = np.clip(expit(th[js] - b[is_]), 1e-12, 1 - 1e-12)
        return -np.sum(ss * np.log(p) + (ns - ss) * np.log(1 - p)) + reg * np.sum(params**2)

    res = minimize(nll, np.zeros(J + F), method='L-BFGS-B', options={'maxiter': 5000})
    theta_l, w = res.x[:J], res.x[J:J + F]
    return theta_l, w, Xn @ w


def build_card_features(card_path, antibiotics):
    """Extract 8 molecular features per antibiotic from CARD aro_categories_index.tsv.

    Features: [n_resistance_genes, n_gene_families,
               efflux, inactivation, target_alteration,
               target_protection, target_replacement, reduced_permeability]

    Drug class mapping uses CARD's own terminology (e.g. 'penicillin' not 'penam',
    'peptide antibiotic' not 'lipopeptide'). All 42 ATLAS antibiotics map to
    non-zero feature vectors. Antibiotics within the same drug class share identical
    feature vectors (15 distinct vectors for 42 antibiotics); intra-class potency
    variation is captured by the free IRT parameters, not the LLTM decomposition.
    """
    aro = pd.read_csv(card_path, sep='\t')
    abx_map = {
        # Aminoglycosides
        'amikacin': ['aminoglycoside'], 'gentamicin': ['aminoglycoside'],
        'tobramycin': ['aminoglycoside'],
        # Penicillins — CARD uses 'penicillin beta-lactam', matched by 'penicillin'
        'ampicillin': ['penicillin'], 'ampicillin sulbactam': ['penicillin'],
        'amoxycillin clavulanate': ['penicillin'], 'oxacillin': ['penicillin'],
        'penicillin': ['penicillin'], 'piperacillin tazobactam': ['penicillin'],
        # Cephalosporins (CARD does not distinguish cephamycins)
        'cefepime': ['cephalosporin'], 'cefiderocol': ['cephalosporin'],
        'cefixime': ['cephalosporin'], 'cefoxitin': ['cephalosporin'],
        'cefpodoxime': ['cephalosporin'], 'ceftaroline': ['cephalosporin'],
        'ceftazidime': ['cephalosporin'], 'ceftazidime avibactam': ['cephalosporin'],
        'ceftibuten': ['cephalosporin'], 'ceftolozane tazobactam': ['cephalosporin'],
        'ceftriaxone': ['cephalosporin'],
        # Carbapenems
        'doripenem': ['carbapenem'], 'ertapenem': ['carbapenem'],
        'imipenem': ['carbapenem'], 'meropenem': ['carbapenem'],
        'meropenem vaborbactam': ['carbapenem'], 'imipenem relebactam': ['carbapenem'],
        # Monobactam
        'aztreonam': ['monobactam'],
        # Fluoroquinolones
        'ciprofloxacin': ['fluoroquinolone'], 'levofloxacin': ['fluoroquinolone'],
        'moxifloxacin': ['fluoroquinolone'],
        # Tetracyclines
        'minocycline': ['tetracycline'], 'tetracycline': ['tetracycline'],
        'tigecycline': ['glycylcycline', 'tetracycline'],
        # Other classes
        'azithromycin': ['macrolide'], 'vancomycin': ['glycopeptide'],
        'linezolid': ['oxazolidinone'],
        'daptomycin': ['peptide antibiotic'],  # CARD has no 'lipopeptide' class
        'colistin': ['peptide antibiotic'],
        'metronidazole': ['nitroimidazole'],
        'trimethoprim sulfa': ['sulfonamide', 'diaminopyrimidine'],
        'nitrofurantoin': ['nitrofuran'], 'chloramphenicol': ['phenicol'],
    }
    mechs = ['antibiotic efflux', 'antibiotic inactivation', 'antibiotic target alteration',
             'antibiotic target protection', 'antibiotic target replacement',
             'reduced permeability to antibiotic']
    X = np.zeros((len(antibiotics), 8))
    for i, abx in enumerate(antibiotics):
        cc = abx_map.get(abx.lower(), [])
        sub = aro[aro['Drug Class'].str.contains('|'.join(cc), case=False, na=False)] if cc else aro.iloc[0:0]
        X[i, 0] = len(sub)
        X[i, 1] = sub['AMR Gene Family'].nunique() if len(sub) else 0
        for k, m in enumerate(mechs):
            X[i, 2 + k] = sub['Resistance Mechanism'].str.contains(m, case=False, na=False).sum() if len(sub) else 0
    return X


def main():
    ap = argparse.ArgumentParser(description="AMR-IRT pipeline (Vivli 2026)")
    ap.add_argument('--atlas', required=True, help='Path to ATLAS CSV extract (from Vivli)')
    ap.add_argument('--card', required=True, help='Path to CARD aro_categories_index.tsv')
    ap.add_argument('--out', default='outputs', help='Output directory')
    args = ap.parse_args()

    if not os.path.exists(args.atlas):
        sys.exit(f"ERROR: ATLAS file not found: {args.atlas}\n"
                 "Obtain it via the Vivli AMR Register (see REPRODUCE.md). "
                 "Raw ATLAS data cannot be redistributed.")
    if not os.path.exists(args.card):
        sys.exit(f"ERROR: CARD file not found: {args.card}\n"
                 "Download from https://card.mcmaster.ca/download")

    os.makedirs(args.out, exist_ok=True)

    print("Loading ATLAS matrix...")
    pathogens, antibiotics, R, N, mask = load_atlas_matrix(args.atlas)
    print(f"  {len(pathogens)} species x {len(antibiotics)} antibiotics, "
          f"{int(mask.sum())}/{mask.size} cells observed ({mask.mean():.1%})")

    print("Fitting 1PL IRT (Rasch)...")
    theta, b, _ = fit_irt_1pl(R, mask, N)

    print("Building CARD features and fitting LLTM...")
    X = build_card_features(args.card, antibiotics)
    theta_l, w, b_l = fit_lltm(R, mask, N, X)
    rho_t = spearmanr(theta, theta_l).correlation
    rho_b = spearmanr(b, b_l).correlation
    print(f"  IRT theta vs LLTM theta: rho = {rho_t:.3f}")
    print(f"  IRT b     vs LLTM b:     rho = {rho_b:.3f}")

    order = np.argsort(-theta)
    avg = np.array([np.nanmean(R[j, mask[j]]) if mask[j].sum() else 0 for j in range(len(pathogens))])
    pd.DataFrame([{'rank': r + 1, 'species': pathogens[j], 'theta': theta[j],
                   'avg_resistance': avg[j], 'n_antibiotics': int(mask[j].sum())}
                  for r, j in enumerate(order)]).to_csv(
        os.path.join(args.out, 'superbug_leaderboard.csv'), index=False)

    bo = np.argsort(-b)
    pd.DataFrame([{'rank': r + 1, 'antibiotic': antibiotics[i], 'b_irt': b[i], 'b_lltm': b_l[i],
                   'n_species': int(mask[:, i].sum())}
                  for r, i in enumerate(bo) if mask[:, i].sum()]).to_csv(
        os.path.join(args.out, 'antibiotic_potency.csv'), index=False)

    print(f"Done. Outputs written to {args.out}/")


if __name__ == '__main__':
    main()
