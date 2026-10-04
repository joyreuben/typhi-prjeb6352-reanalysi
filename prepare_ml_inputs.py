#!/usr/bin/env python3
"""Build the Random Forest inputs for the 157 core Typhi samples.

SNP matrix: from the whole-genome alignment (gubbins/clean.full.aln), with each sample's
Gubbins recombination blocks set to N, keeping biallelic sites that are called (A/C/G/T)
in all 157 samples and variable among them. 1 = minor allele, 0 = major allele.
Columns are CT18 positions, so position-based exclusions (e.g. an LD window) stay possible.

Labels: Mykrobe calls. MDR = resistant to ampicillin + chloramphenicol + co-trimoxazole.

Usage (typhi_amr env):  python prepare_ml_inputs.py
Output: ml/snp_matrix.tsv.gz, ml/sites.tsv, ml/labels.tsv
"""
import os
import re

import numpy as np
import pandas as pd

ALN = "gubbins/clean.full.aln"
GFF = "gubbins/gubbins.recombination_predictions.gff"
MYKROBE = "mykrobe/mykrobe_out_predictResults.tsv"
SAMPLES = "results/typhi_core_runs.txt"
CHUNK = 500_000

os.makedirs("ml", exist_ok=True)
samples = open(SAMPLES).read().split()

# --- read alignment (one uint8 row per sample) ---
seqs, name, buf = {}, None, []
with open(ALN) as fh:
    for line in fh:
        if line.startswith(">"):
            if name in samples:
                seqs[name] = np.frombuffer("".join(buf).upper().encode(), dtype=np.uint8).copy()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name in samples:
        seqs[name] = np.frombuffer("".join(buf).upper().encode(), dtype=np.uint8).copy()
missing = set(samples) - set(seqs)
assert not missing, f"not in alignment: {missing}"
aln = np.vstack([seqs[s] for s in samples])
del seqs
L = aln.shape[1]
print(f"alignment: {aln.shape[0]} samples x {L:,} bp")

# --- mask each sample's recombination blocks ---
idx = {s: i for i, s in enumerate(samples)}
n_blocks = masked = 0
for line in open(GFF):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    start, end = int(f[3]), int(f[4])
    taxa = re.search(r'taxa="([^"]*)"', f[8]).group(1).split()
    rows = [idx[t] for t in taxa if t in idx]
    if rows:
        aln[rows, start - 1:end] = ord("N")
        n_blocks += 1
        masked += len(rows) * (end - start + 1)
print(f"recombination: masked {n_blocks} blocks ({masked:,} sample-bp)")

# --- find core biallelic variable sites, chunk by chunk to keep memory modest ---
bases = np.array([ord(b) for b in "ACGT"], dtype=np.uint8)
pos_list, major_list, minor_list, cols = [], [], [], []
for a in range(0, L, CHUNK):
    blk = aln[:, a:a + CHUNK]
    counts = np.stack([(blk == b).sum(axis=0) for b in bases])      # 4 x chunk
    core = counts.sum(axis=0) == blk.shape[0]                         # no N/- in any sample
    biallelic = (counts > 0).sum(axis=0) == 2
    for j in np.where(core & biallelic)[0]:
        c = counts[:, j]
        order = np.argsort(c)[::-1]
        maj, mnr = bases[order[0]], bases[order[1]]
        pos_list.append(a + j + 1)
        major_list.append(chr(maj))
        minor_list.append(chr(mnr))
        cols.append((blk[:, j] == mnr).astype(np.uint8))

X = pd.DataFrame(np.column_stack(cols), index=samples, columns=pos_list)
X.index.name = "genome"
sites = pd.DataFrame({"pos": pos_list, "major": major_list, "minor": minor_list,
                      "minor_count": X.sum(axis=0).values})
print(f"core biallelic SNP sites among {len(samples)} samples: {X.shape[1]:,}")
print(f"  singletons (minor allele in 1 sample): {(sites.minor_count == 1).sum():,}")

# --- labels from Mykrobe ---
m = pd.read_csv(MYKROBE, sep="\t").set_index("genome").loc[samples]
R = lambda col: m[col].astype(str).str.startswith("R")
lab = pd.DataFrame(index=pd.Index(samples, name="genome"))
lab["genotype"] = m["genotype"]
lab["MDR"] = (R("ampicillin") & R("chloramphenicol") & R("trimethoprim-sulfamethoxazole")).astype(int)
lab["IncHI1"] = (m["IncHI1A"] > 0).astype(int)
lab["gyrA_S83"] = ((m["gyrA_S83F"] + m["gyrA_S83Y"]) > 0).astype(int)
lab["batch"] = lab.index.str[:6]
lab["read_len"] = np.where(lab.batch == "ERR573", 150, 100)
print(f"MDR: {lab.MDR.sum()} / {len(lab)}")

X.to_csv("ml/snp_matrix.tsv.gz", sep="\t", compression="gzip")
sites.to_csv("ml/sites.tsv", sep="\t", index=False)
lab.to_csv("ml/labels.tsv", sep="\t")
print("wrote ml/snp_matrix.tsv.gz, ml/sites.tsv, ml/labels.tsv")
