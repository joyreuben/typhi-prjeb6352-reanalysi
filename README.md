# Nigerian *Salmonella* Typhi: genomic epidemiology and leakage-aware AMR prediction

A corrected re-analysis of the study *"Genomic Epidemiology and Machine Learning Classification of
Fluoroquinolone Resistance-Associated Sublineages in Salmonella Typhi from Nigeria"*, using all
229 Nigerian sequencing runs in ENA project [PRJEB6352](https://www.ebi.ac.uk/ena/browser/view/PRJEB6352)
(Wellcome Sanger Institute, Illumina paired-end).

The original study trained a Random Forest on genome-wide SNPs to predict *parC* S80R and reported
AUROC 1.00 with random cross-validation. In a clonal pathogen, random CV lets lineage-linked SNPs
stand in for the label (population-structure leakage). This re-analysis adds recombination
filtering, phylogeny-aware cross-validation and label-permutation tests.

**Status:** genomics and the first ML analysis are complete; write-up in progress.

## Findings so far

1. **70 of the 229 runs are not Typhi**, although ENA labels all of them *S.* Typhi.
   They have ~48,000 SNPs against the CT18 reference (Typhi runs: 312–861), cover only ~87% of
   CT18, and lack the Typhi-specific *tviA/tviB* (Vi capsule) and *staA* genes. SeqSero2 assigns
   them to other serovars: Typhimurium 29, Enteritidis 17, Dublin 9, and 11 others
   (`results/seqsero2_all.tsv`).
2. **2 Typhi runs look mixed** (ERR731385, ERR657450: 79,972 and 27,374 heterozygous calls vs a
   median of ~630) and were excluded. **Final set: 157 Typhi genomes.**
3. **No *parC* S80R in any of the 229 runs** (Typhi Mykrobe and Snippy agree). The only quinolone
   resistance determinants are *gyrA* S83Y (6), *gyrA* S83F (1) and *qnrS1* (1).
4. **The isolates are not a single clone.** GenoTyphi genotypes: 3.1.1 (98), 2.3.1 (17), 2.2 (15),
   0.0.3 (8), 4.1 (8), 2.1 (4), 3.3 (2), 2.3.2 (2), 0.0.1 (2), 4.1.1 (1).
5. **MDR is lineage-bound.** 74/157 are multidrug resistant (ampicillin, chloramphenicol,
   co-trimoxazole: *blaTEM-1D*, *catA1*, *dfrA*/*sul*). All 74 are genotype 3.1.1 and carry an IncHI1
   plasmid. No ceftriaxone, azithromycin or carbapenem resistance was found.
6. **Low recombination:** Gubbins found 25 recombination blocks (1,027 SNPs) vs 3,185 SNPs by mutation.

Because *parC* S80R is absent, the ML target is **MDR status**. MDR genes are plasmid-borne, so any
chromosomal SNP signal reflects lineage, which makes it a clean test case for leakage.

## Pipeline

| Step | Tool | Script / command | Output |
|---|---|---|---|
| 1. Download | wget (ENA); SRA toolkit for 2 runs with broken ENA files | `download_reads.py` | `raw_reads/` (not in repo) |
| 2. QC + trimming | FastQC, MultiQC, Trimmomatic 0.40 | `trim_all.sh` | `qc/multiqc_raw/multiqc_raw.html`, `trimmed/trim_summary.tsv` |
| 3. Variant calling | Snippy 4.6.0 vs CT18 (AL513382.1, chromosome) | `snippy_all.sh` | `results/snippy_summary.tsv` |
| 4. Species check | SeqSero2 (k-mer mode) | see below | `results/seqsero2_all.tsv`, `results/typhi_runs.txt` |
| 5. Core alignment | snippy-core `--mask auto` (157 Typhi) | see below | `core/` |
| 6. Recombination | Gubbins 3.3 | see below | `gubbins/` |
| 7. Phylogeny | IQ-TREE (ModelFinder + ASC, 1000 UFBoot) | see below | `tree/typhi.treefile` |
| 8. Genotype, QRDR, AMR, plasmids | Typhi Mykrobe (panel 20240407) + GenoTyphi parser | `mykrobe_all.sh` | `mykrobe/mykrobe_out_predictResults.tsv` |
| 9. ML inputs | Python | `prepare_ml_inputs.py` | `ml/snp_matrix.tsv.gz`, `ml/labels.tsv` |
| 10. Random Forest | scikit-learn | `ml_mdr.py` | `ml/cv_results.tsv` |

Trimmomatic settings: `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50`.

Commands for the steps without a script:

```bash
# 4. serovar (env with seqsero2)
SeqSero2_package.py -m k -t 2 -n $r -d results/seqsero2/$r -i trimmed/$r/${r}_1P.fastq.gz trimmed/$r/${r}_2P.fastq.gz

# 5. core alignment of the 157 Typhi
snippy-core --ref ref/CT18.gbk --mask auto --prefix core/core $(sed 's#^#snippy/#' results/typhi_core_runs.txt)

# 6. recombination
snippy-clean_full_aln core/core.full.aln > gubbins/clean.full.aln
cd gubbins && run_gubbins.py --threads 4 --prefix gubbins clean.full.aln

# 7. tree from recombination-filtered core SNPs
snp-sites -c -o tree/core_snps.aln gubbins/gubbins.filtered_polymorphic_sites.fasta
iqtree -s tree/core_snps.aln -m MFP+ASC -B 1000 -T 4 --prefix tree/typhi
```

### ML design (`ml_mdr.py`)

The same Random Forest (500 trees, balanced class weights) on core SNPs found in at least 2 samples,
evaluated four ways:

- **random:** stratified 5-fold CV, repeated 5 times (the original study's approach)
- **clade:** leave-clade-out, with 10 clades cut from the IQ-TREE tree by patristic distance
- **genotype:** leave-genotype-out
- **baseline:** "is the isolate genotype 3.1.1?" as the only predictor

Each CV scheme gets a label-permutation p-value. The whole analysis is repeated within genotype 3.1.1
only, asking whether the chromosome predicts which 3.1.1 isolates acquired the MDR plasmid.

## ML results (`ml/cv_results.tsv`)

Random Forest predicting MDR from chromosomal core SNPs; AUROC on pooled out-of-fold predictions;
permutation p from 100 label shuffles (0.010 is the smallest possible value).

| Samples | Evaluation | AUROC | Balanced acc. | Fraction scored | Perm. p |
|---|---|---|---|---|---|
| all 157 | baseline rule: genotype == 3.1.1 | 0.855 | 0.855 | 1.00 | – |
| all 157 | random 5-fold CV | 0.924 | 0.913 | 1.00 | 0.010 |
| all 157 | leave-clade-out | not estimable | – | 0.24 | – |
| all 157 | leave-genotype-out | not estimable | – | 0.38 | – |
| 98 genotype 3.1.1 | random 5-fold CV | 0.771 | 0.809 | 1.00 | 0.010 |
| 98 genotype 3.1.1 | leave-clade-out (10 Ward clades) | 0.707 | 0.599 | 1.00 | 0.010 |

- **Across all isolates, the model mostly learns lineage.** Random CV (0.924) is only a little above
  the one-line rule "is it genotype 3.1.1?" (0.855). The most important SNPs are alleles found only
  in non-MDR isolates, i.e. markers of the non-MDR lineages.
- **Lineage-held-out performance cannot be measured across all isolates.** Every MDR isolate is in
  genotype 3.1.1, so when that lineage is held out the training data has no MDR examples. A claim of
  "AUROC 1.00" for this kind of data only shows that the model recognises lineages.
- **Within genotype 3.1.1, a modest signal survives holding out sub-clades** (AUROC 0.707, above
  all 100 permutations, max null 0.678), but balanced accuracy drops to 0.599. This fits the IncHI1
  MDR plasmid being inherited along 3.1.1 sub-branches: chromosomal SNPs mark the sub-lineages that
  carry it. It is still population structure at a finer scale, not a resistance mechanism.

## What is not in this repo

Reads (~115 GB raw + trimmed), per-sample Snippy folders, whole-genome alignments (`core.full.aln`,
`gubbins/clean.full.aln`, ~750 MB each) and the CT18 GenBank file. All of these can be regenerated:
reads from PRJEB6352 using `prjeb6352_final_manifest.csv` (checksums in `ena_md5_all.txt`), and CT18
from NCBI accession AL513382.1.

`CLAUDE.md` is the dated lab notebook for the project, including problems hit along the way and how
they were resolved.

## Environment

WSL2 Ubuntu 24.04, conda. Main env `typhi_amr` (snippy, trimmomatic, fastqc, gubbins, iqtree, snp-sites,
python 3.10, pandas, scikit-learn, dendropy). Separate envs: `seqsero2`, `mykrobe`, `multiqc`.
Gubbins 3.3 and Mykrobe need `setuptools<81` (they import `pkg_resources`).
