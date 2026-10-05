#!/usr/bin/env python3
"""Random Forest: can chromosomal core SNPs predict MDR, and how much of that is just lineage?

Compares cross-validation schemes on the same model and data:
  random  - stratified 5-fold, repeated (what the original study did)
  clade   - leave-clade-out: GroupKFold on clusters cut from the IQ-TREE tree (patristic distance,
            Ward linkage; average/complete linkage gave one giant cluster + singletons)
  genotype- leave-genotype-out: GroupKFold on GenoTyphi genotype
Plus:
  baseline    - "is it genotype 3.1.1?" as the only predictor (pure lineage)
  permutation - labels shuffled N times, same CV, to give an empirical p-value
  within 3.1.1- the same comparison restricted to the 98 genotype-3.1.1 samples
                (does the chromosome predict who acquired the MDR plasmid inside one lineage?)

AUROC is computed on pooled out-of-fold predictions. If a training fold lacks one class (e.g. all MDR
isolates are in the held-out clade), its test samples cannot be scored: they are left out and the
fraction scored is reported. AUROC is NaN ("not estimable") when the scored samples hold one class.

Usage (typhi_amr env):  python ml_mdr.py [n_permutations]   (default 100)
Output: ml/cv_results.tsv, ml/permutation_null.tsv, ml/feature_importance.tsv, ml/clades.tsv
"""
import sys

import dendropy
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold

N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 100
N_CLADES = 10          # clusters cut from the tree for leave-clade-out CV
N_REPEATS = 5          # repeats of random 5-fold CV
SEED = 1
RF = dict(n_estimators=500, class_weight="balanced", n_jobs=4, random_state=SEED)
RF_PERM = dict(RF, n_estimators=200)

X_all = pd.read_csv("ml/snp_matrix.tsv.gz", sep="\t", index_col=0)
lab = pd.read_csv("ml/labels.tsv", sep="\t", index_col=0).loc[X_all.index]
X_all = X_all.loc[:, X_all.sum(axis=0) >= 2]          # drop singletons
print(f"{X_all.shape[0]} samples, {X_all.shape[1]} SNPs (minor allele in >=2 samples), "
      f"MDR {lab.MDR.sum()}/{len(lab)}")

# --- clades from the IQ-TREE tree (patristic distance, average linkage) ---
tree = dendropy.Tree.get(path="tree/typhi.treefile", schema="newick", preserve_underscores=True)
pdm = tree.phylogenetic_distance_matrix()
taxa = {t.label: t for t in tree.taxon_namespace}
names = list(X_all.index)
D = np.array([[pdm.distance(taxa[a], taxa[b]) if a != b else 0.0 for b in names] for a in names])


def tree_clades(sub, k):
    idx = [names.index(s) for s in sub]
    Z = linkage(squareform(D[np.ix_(idx, idx)], checks=False), method="ward")
    return pd.Series(fcluster(Z, t=k, criterion="maxclust"), index=sub)


lab["clade"] = tree_clades(names, N_CLADES)
lab[["genotype", "clade", "MDR"]].to_csv("ml/clades.tsv", sep="\t")
print("\nclade x MDR:")
print(pd.crosstab([lab.clade], lab.MDR).rename(columns={0: "non-MDR", 1: "MDR"}).to_string())


def splits(scheme, y, groups, seed):
    if scheme == "random":
        return list(StratifiedKFold(5, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))
    n = min(5, len(np.unique(groups)))
    return list(GroupKFold(n).split(np.zeros(len(y)), y, groups))


def oof_scores(X, y, folds, rf):
    p = np.full(len(y), np.nan)
    for tr, te in folds:
        if len(np.unique(y[tr])) < 2:          # training fold has one class only: cannot score
            continue
        m = RandomForestClassifier(**rf).fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p


def evaluate(X, y, groups_by_scheme, rf, n_rep):
    out = {}
    for scheme, groups in groups_by_scheme.items():
        reps = n_rep if scheme == "random" else 1
        aucs, bals, fracs = [], [], []
        for r in range(reps):
            p = oof_scores(X, y, splits(scheme, y, groups, SEED + r), rf)
            ok = ~np.isnan(p)
            fracs.append(ok.mean())
            if len(np.unique(y[ok])) < 2:
                aucs.append(np.nan); bals.append(np.nan)
                continue
            aucs.append(roc_auc_score(y[ok], p[ok]))
            bals.append(balanced_accuracy_score(y[ok], p[ok] >= 0.5))
        out[scheme] = (np.mean(aucs), np.mean(bals), np.mean(fracs))
    return out


rows, null_rows = [], []
rng = np.random.default_rng(SEED)
for subset_name, keep in [("all 157", lab.index), ("within 3.1.1", lab.index[lab.genotype == "3.1.1"])]:
    Xs = X_all.loc[keep]
    Xs = Xs.loc[:, (Xs.sum(axis=0) >= 2) & (Xs.sum(axis=0) <= len(Xs) - 2)].values
    y = lab.loc[keep, "MDR"].values
    sub_clades = tree_clades(list(keep), N_CLADES) if subset_name != "all 157" else lab.loc[keep, "clade"]
    if subset_name != "all 157":
        print(f"\nclade x MDR ({subset_name}):")
        print(pd.crosstab(sub_clades, y).rename(columns={0: "non-MDR", 1: "MDR"}).to_string())
    schemes = {"random": None, "clade": sub_clades.values}
    if subset_name == "all 157":
        schemes["genotype"] = lab.loc[keep, "genotype"].values
    print(f"\n== {subset_name}: {len(y)} samples, {Xs.shape[1]} SNPs, MDR {y.sum()}")

    real = evaluate(Xs, y, schemes, RF, N_REPEATS)
    if subset_name == "all 157":
        base = (lab.loc[keep, "genotype"] == "3.1.1").astype(int).values
        rows.append(dict(subset=subset_name, scheme="baseline: genotype==3.1.1",
                         AUROC=roc_auc_score(y, base), bal_acc=balanced_accuracy_score(y, base),
                         frac_scored=1.0, perm_p=np.nan))

    null = {s: [] for s in schemes}
    for i in range(N_PERM):
        yp = rng.permutation(y)
        res = evaluate(Xs, yp, schemes, RF_PERM, 1)
        for s in schemes:
            null[s].append(res[s][0])
            null_rows.append(dict(subset=subset_name, scheme=s, perm=i, AUROC=res[s][0]))
        if (i + 1) % 10 == 0:
            print(f"  permutations {i + 1}/{N_PERM}", flush=True)

    for s, (auc, bal, frac) in real.items():
        nl = [a for a in null[s] if not np.isnan(a)]
        p = np.nan if np.isnan(auc) or not nl else (1 + sum(a >= auc for a in nl)) / (1 + len(nl))
        rows.append(dict(subset=subset_name, scheme=s, AUROC=auc, bal_acc=bal, frac_scored=frac, perm_p=p))

res = pd.DataFrame(rows)
res.to_csv("ml/cv_results.tsv", sep="\t", index=False, float_format="%.3f")
pd.DataFrame(null_rows).to_csv("ml/permutation_null.tsv", sep="\t", index=False, float_format="%.4f")

# feature importance from a model fit on all 157 (descriptive only)
X = X_all.values
m = RandomForestClassifier(**RF).fit(X, lab.MDR.values)
fi = pd.DataFrame({"pos": X_all.columns.astype(int), "importance": m.feature_importances_,
                   "minor_count": X.sum(axis=0),
                   "minor_in_MDR": X[lab.MDR.values == 1].sum(axis=0)})
fi.sort_values("importance", ascending=False).to_csv("ml/feature_importance.tsv", sep="\t", index=False)

print("\n" + res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\nwrote ml/cv_results.tsv, ml/permutation_null.tsv, ml/feature_importance.tsv, ml/clades.tsv")
