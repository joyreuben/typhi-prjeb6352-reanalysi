# Typhi study — Nigerian *S.* Typhi fluoroquinolone resistance (redo)

## What this project is
A corrected redo of the study "Genomic Epidemiology and Machine Learning Classification of
Fluoroquinolone Resistance-Associated Sublineages in *Salmonella* Typhi from Nigeria".

The original used 104 isolates from ENA **PRJEB6352**, Snippy vs CT18, ABRicate/CARD, genotyphi,
and a Random Forest on 2,265 SNPs to predict **parC S80R** status (parC locus excluded,
5-fold stratified CV). It reported AUROC 1.00, which is **population-structure leakage**: all
isolates are one clone (genotype 3.2.1), parC carriers form their own sub-clade (they separate
fully on PC1), so linked SNPs across the genome still encode the label.

This redo uses **all 229 Nigerian runs** in PRJEB6352 (Wellcome Sanger, 2008–2013, paired-end, ~62 GB).

## Pipeline
1. Download reads — `download_reads.py` → `raw_reads/<run>/` (log: `download.log`)
2. QC — FastQC + Trimmomatic
3. Variant calling — Snippy vs CT18 reference
4. Recombination filtering — Gubbins *(new vs original)*
5. Phylogeny — IQ-TREE
6. Genotyping — genotyphi (`genotyphi/`)
7. AMR genes — ABRicate with CARD
8. QRDR mutations (gyrA, gyrB, parC, parE) extracted from VCFs
9. SNP matrix with an **LD-window exclusion around parC**, not just the locus *(new)*
10. Random Forest with **phylogeny/clade-aware CV** and a **label-permutation test** *(new)*;
    report naive random CV and clade-aware CV side by side so the gap shows the leakage.

## Environment
- WSL Ubuntu 24.04, 8 CPU, 7.6 GB RAM (keep memory-heavy steps modest)
- Conda env: `conda activate typhi_amr` (snippy, trimmomatic, fastqc, abricate, gubbins, iqtree, seqkit, python 3.10, pandas, scikit-learn)
- CT18 reference downloaded (`ref/CT18.gbk`). Not yet confirmed: `abricate --setupdb` run
- **New conda envs get setuptools 84 (no `pkg_resources`)** → older tools (gubbins, mykrobe) crash; fix with
  `conda install -c conda-forge -c bioconda "setuptools<81" -y` in that env. Other envs: `seqsero2`, `mykrobe`, `multiqc`.

## Files
- `prjeb6352_final_manifest.csv` — 229 runs: run_accession, sample_accession, fastq_ftp
- `download_reads.py` — parallel (8) wget downloader; skips finished files, resumes `.part` files
- `raw_reads/` — DELETED 2026-10-02 to save space (FastQC + trimming done; `trimmed/` is the input now)
- `genotyphi/` — clone of github.com/katholt/genotyphi

## Working style
- The user runs pipeline commands themself. Give copyable commands with a simple explanation,
  one step at a time, and wait for the pasted output. Read-only status checks are fine.
- Explain in simple terms; don't do everything at once.
- A PC restart kills background jobs in WSL; relaunch with
  `nohup python3 -u download_reads.py > download.log 2>&1 &` (safe to rerun).

## Status log
- 2026-09-18: project planned, env created, download started
- 2026-09-24: download at 451/458, last 7 files in progress.
- 2026-09-25: download finished, but `gzip -t` found 30 corrupt files (20 runs); all failed ENA md5.
  Deleted and re-downloading (`download2.log`). Checksums: `ena_md5_all.txt`, `redownload_md5.txt`.
  **ENA server bug:** ERR573956_2 and ERR573940_2 are empty directories on ENA (HTTP returns an
  HTML listing), so both mates of those 2 runs are fetched from NCBI SRA via prefetch/fasterq-dump
  (amr env). Verify those by read count, not ENA md5: ERR573940 = 1,337,552 pairs (2,675,104 reads),
  ERR573956 = 1,316,274 pairs (2,632,548 reads). SRA dump in `sra_fetch/fastq/` matched both.
  Re-download done; 28/28 ENA-fetched files pass md5 (`redownload_md5_check.txt`).
  SRA fetch done: both mates of ERR573940 + ERR573956 replaced with gzipped SRA dumps (gzip -t OK).
  Final check: 229/229 manifest runs present, 458 fastq.gz, no mates missing, no files <5 MB, 63 GB.
  `sra_fetch/` (uncompressed FASTQs + .sra) can be deleted. Download stage complete.
- 2026-09-25: FastQC raw done (`qc/fastqc_raw/`, 458/458). MultiQC in its own env `multiqc`
  (typhi_amr can't solve it; create with `-c conda-forge -c bioconda` order — strict priority).
  Report: `qc/multiqc_raw/multiqc_raw.html`. Data is clean: quality/adapters/dup pass on all,
  coverage 38–177x (median 95x). GC "warn" on all = normal for bacteria. Per-tile fail on one lane
  (ERR573962–973 _1). **Read length batch: ERR573xxx = 150 bp; ERR657xxx/ERR731xxx = 100 bp** —
  check as batch-effect confounder vs parC label. **ERR657395 GC 54% (others 50–52%)** — possible
  contamination; check CT18 mapping rate after Snippy.
- 2026-09-27: Trimmomatic done, 229/229, 0 failed (`trim_all.sh`, rerunnable; settings
  ILLUMINACLIP TruSeq3-PE-2 :2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50).
  Output `trimmed/<run>/<run>_{1,2}P.fastq.gz` (use P only), summary `trimmed/trim_summary.tsv`.
  Both-surviving 68–94% (median 89%), dropped ≤5.6%; est. coverage after trim min ~30x, median ~81x.
  Lowest survival are 150 bp runs (ERR573959 68%, ERR573948 71%) — R2 tails trimmed, still ≥36x.
  Lowest coverage ERR657451 ~30x.
- 2026-09-27: CT18 reference `ref/CT18.gbk` = AL513382.1 chromosome only (no plasmids), 4,809,037 bp,
  4,600 CDS (efetch threw curl 56 at the end but file is complete). **In this record `/gene` is the
  STY locus tag; real names are in `/gene_synonym`**, so Snippy output will say STY3351, not parC.
  QRDR genes (all but gyrB on minus strand): gyrA STY2499 complement(2331373..2334009);
  parC STY3351 complement(3194450..3196708); parE STY3359 complement(3201464..3203356);
  gyrB STY3943 3808932..3811346.
- 2026-09-27: Snippy 4.6.0 test ERR573926: 99.05% mapped, 69x mean, 1.3% of genome <10x, 793 variants
  (709 SNP), QRDR genes 70x with no variants. Contig name in Snippy output is `AL513382`.
  samtools is in typhi_amr bin. Full run: `snippy_all.sh` (2 parallel x 4 CPU, rerunnable).
  **Next: run snippy_all.sh, then per-sample mapping/depth summary (check ERR657395), snippy-core.**
- 2026-09-30: snippy_all.sh stopped at 41/229 (restart); relaunched → `snippy2.log`, 75/229 at 20:08, 0 failed.
  **ERR657395 fails: 35% mapped to CT18 (others ~97–99%), 47,505 variants (others ~600–800)** →
  contaminated/mixed, likely exclude from snippy-core (confirm with Kraken2 or a species check if needed).
- 2026-10-02: Snippy stopped at 177/229 done (ERR731387/388 cut mid-run; rerun redoes them). Windows C: had filled
  (WSL disk `C:\Users\Joy\AppData\Local\wsl\{aeb32ab6-...}\ext4.vhdx` 191 GB). On Oct 1 (Copilot chat) BAMs of the
  177 done runs were deleted and key files copied to `results/snippy_keep/` (5.3 GB, duplicate of snippy/).
  Deleted `raw_reads/` (FastQC+trim done; re-downloadable). fstrim + diskpart compact → C: 106 GB free.
  **If C: fills again:** Linux-side deletes don't free C: until fstrim + `wsl --shutdown` + diskpart `compact vdisk`.
  Next: rerun snippy_all.sh for the last 52 runs (trimmed reads present for all 52).
- 2026-10-03: **Snippy done 229/229, 0 failed** (`snippy3.log`). Per-sample summary `results/snippy_summary.tsv`.
  **70 runs are NOT Typhi** despite ENA labelling all 229 as Typhi (`results/ena_taxa.tsv`): ~48k variants vs CT18
  (Typhi runs: 312–~1,000, median ~600), only ~87% of CT18 covered, and they lack the Typhi-specific genes
  tviA/tviB (Vi capsule, SPI-7) and staA. Looks like another Salmonella serovar (to confirm with SeqSero2/SISTR).
  Lists: `results/non_typhi_runs.txt` (70, incl. ERR657395; ERR657386 only 18% covered) and `results/typhi_runs.txt` (159).
  Plan: snippy-core etc. on the 159 Typhi only.
- 2026-10-03: SeqSero2 (env `seqsero2`, k-mer mode, `results/seqsero2/`) confirms: ERR573926 = Typhi (control);
  ERR573928 = Agama, ERR657376 = Welikade, ERR657395 = Typhimurium, ERR731377 = Durban. So the 70 are a MIX of
  non-Typhi serovars (similar ~48k SNP counts only because any non-Typhi serovar is about equally far from CT18).
- 2026-10-03: snippy-core (`--mask auto`, 32 repeat regions/202 kb masked) → `core/`. Excluded 2 likely mixed
  Typhi samples with huge het counts: ERR731385 (79,972 het) and ERR657450 (27,374) vs median ~630.
  Final set **157 Typhi** (`results/typhi_core_runs.txt`), 2,063 core SNP sites (+ CT18 Reference row).
  SeqSero2 on all 70 non-Typhi running in background → `results/seqsero2_all.tsv`.
  Gubbins 3.3.4 crashes: `No module named 'pkg_resources'` (setuptools 84 removed it) → fix with setuptools<81.
- 2026-10-03: Gubbins done (8.5 min, `gubbins/`, 158 seqs incl. Reference; hit max 5 iterations). 25 recombination
  blocks, 1,027 SNPs in recombination vs 3,185 outside (low, as expected for Typhi); biggest on Reference branch
  (4 blocks, 275 SNPs) and Node_140. SeqSero2 all 70 done (`results/seqsero2_all.tsv`), none Typhi: Typhimurium 29,
  Enteritidis 17, Dublin 9, Paratyphi C/Choleraesuis 3, I 4,[5],12:i:- 3, Virchow 2, Agama 2, 1 each Welikade,
  Saintpaul, Durban, I 35:y:l,w, and 1 untypeable. Next: IQ-TREE on Gubbins-filtered core SNPs.
- 2026-10-03: IQ-TREE done (`tree/typhi.treefile`, prefix `tree/typhi`): input = Gubbins-filtered core SNPs
  (`tree/core_snps.aln`, snp-sites -c), 158 seqs, 1,861 SNP sites (809 parsimony-informative). Best model
  K3P+ASC (BIC), UFBoot 1000: 117/124 internal nodes ≥95, all ≥70. Tree is unrooted (incl. CT18 Reference).
  BAMs were deleted for 177 runs, so genotyping plan = Typhi Mykrobe on trimmed reads (genotyphi README recommends it).
- 2026-10-04: Typhi Mykrobe done 159/159 (`mykrobe_all.sh`, env `mykrobe`, JSONs in `mykrobe/`, table
  `mykrobe/mykrobe_out_predictResults.tsv`; ERR573932 failed 1st pass from a parallel skeleton-build race, OK on rerun;
  keep `mykrobe/data/`). All 159 serovar typhi, 157 strong confidence (weak: ERR731385 [mixed], ERR657385).
  **Genotypes (157 core) are NOT one clone 3.2.1:** 3.1.1 = 98, 2.3.1 = 17, 2.2 = 15, 0.0.3 = 8, 4.1 = 8, 2.1 = 4,
  3.3 = 2, 2.3.2 = 2, 0.0.1 = 2, 4.1.1 = 1.
  **parC S80R = 0/229** (Mykrobe and Snippy agree; no parC codon-80 change in Typhi or the 70 non-Typhi either).
  Only QRDR changes: gyrA S83Y 6 + S83F 1 (7/157, all 3.1.1, scattered across batches); qnrS1 1 (ERR731404).
  Snippy also sees gyrA D538N (9), V328I (2), P864L (1), parE A364V (2) — not QRDR-known.
  MDR (amp+chl+sxt; blaTEM-1D, catA1, dfrA/sul) = 74/157, ALL genotype 3.1.1, all IncHI1 plasmid (+18 IncHI1 non-MDR).
  No ceftriaxone/azithro/carbapenem resistance. **The original study's label (parC S80R) does not exist in this data**
  → ML target must be redefined (decision pending with user).
- 2026-10-04: **Decision (user): ML target = MDR status** (74 MDR / 83 non-MDR of 157). Expect naive CV ≈ perfect and
  clade-aware CV to collapse, because MDR = 3.1.1 + IncHI1 only. MDR genes are on the plasmid (not in the CT18
  chromosome core SNPs), so any chromosomal signal is pure lineage. ABRicate on Snippy consensus would miss plasmid
  genes (CT18 ref is chromosome-only) → Mykrobe AMR calls are the AMR result.
- 2026-10-04: `prepare_ml_inputs.py` → `ml/` (snp_matrix.tsv.gz 157 x 1,689 core biallelic SNPs after masking 21
  Gubbins blocks, 881 singletons; sites.tsv with CT18 positions; labels.tsv: genotype, MDR, IncHI1, gyrA_S83, batch,
  read_len). Pipeline step 9's parC LD-window no longer applies (target is plasmid-borne MDR, no chromosomal locus).
  `ml_mdr.py [n_perm]`: RF (500 trees, balanced), singletons dropped; random 5-fold x5 vs leave-clade-out (10 tree
  clusters, patristic/average linkage) vs leave-genotype-out; genotype==3.1.1 baseline; label permutations;
  repeated within 3.1.1 only. Outputs ml/cv_results.tsv etc.
- 2026-10-05: ml_mdr.py run 1 (100 perms): baseline genotype==3.1.1 AUROC 0.855; random CV 0.924 (p=0.01);
  within 3.1.1: random 0.771 (p=0.01) vs clade 0.561 (p=0.17). **Run 1 clade/genotype results invalid:** average
  linkage made one giant clade (121/157 all-sample, 79/98 within 3.1.1), and folds whose training set had no MDR were
  scored as constant 0 → AUROC 0.28/0.15 (artifact, below 0.5). Fixed: Ward linkage (3.1.1 clades 48,12,11,7,5,5,4,4,1,1),
  unscorable folds left out (NaN, `frac_scored` column). All-157 leave-clade/genotype-out is inherently not estimable
  (all MDR in one lineage). Top RF features are minor alleles found only in non-MDR isolates = lineage markers.
  Git: local repo initialised, `.gitignore` + `README.md` written, first commit 245f5af (no remote yet).
