# Onco-eQTM: Pan-Cancer Expression Quantitative Trait Methylation Pipeline

> Manuscript: *Expression Quantitative Trait Methylation Across Multiple Cancer Types with Functional and Therapeutic Characterization using Onco-eQTM*
>
> Authors: Bhanu Teja Korra¹, Mayilaadumveettil Nishana², Rahul Kumar¹
>
> ¹ Department of Biotechnology, Indian Institute of Technology Hyderabad, 502285, India
> ² School of Biology, Indian Institute of Science Education and Research Thiruvananthapuram, Kerala 695551, India
>
> Correspondence: rahulk@bt.iith.ac.in

---

## Overview

This repository contains all analysis scripts used to build the Onco-eQTM database — a pan-cancer resource that maps cis-regulatory DNA methylation to gene expression, drug sensitivity, immune infiltration, miRNA regulation, and pathway activity across 27 TCGA cancer types and 6,880 primary tumor samples.

The pipeline identifies high-confidence CpG–gene associations (eQTMs) using a consensus strategy combining Torch-eCpG and MatrixEQTL, then integrates results with four additional functional layers:

- cis-eQTM associations (CpG ↔ gene expression) — 5.25 million
- Drug response associations (251 compounds, GDSC) — 14.45 million
- Pathway activity associations (1,387 PARADIGM pathways) — 13.6 million
- Immune infiltration associations (68 immune cell types) — 3.55 million
- miRNA-associated CpG associations — 4.52 million

The web database is freely available at: https://project.iith.ac.in/cgntlab/OncoeQTM/

---

## Citation

This manuscript is currently under review. If you use these scripts or data, please cite the Zenodo deposit:

> Korra BT, Nishana M, Kumar R. Expression Quantitative Trait Methylation Across Multiple Cancer Types with Functional and Therapeutic Characterization using Onco-eQTM. [Journal, Year, DOI — to be updated upon publication]

```bibtex
@software{korra_2026_onco_eqtm,
  author    = {Korra, Bhanu Teja and Nishana, Mayilaadumveettil and Kumar, Rahul},
  title     = {Onco-eQTM: Scripts for a Pan-Cancer Atlas of CpG
               Methylation--Gene Expression Associations},
  year      = {2026},
  publisher = {Zenodo},
  version   = {1.0.0},
  doi       = {10.5281/zenodo.XXXXXXX},
  url       = {https://doi.org/10.5281/zenodo.XXXXXXX},
  note      = {Manuscript under review}
}
```

Replace XXXXXXX with your Zenodo record ID after deposit.

---

## Repository Contents

### Download Scripts

These Bash scripts fetch all input datasets from public repositories. They are located in scripts/download/ and can be run individually or all at once using scripts/run_downloads.sh.

- 01_download_clinical_data.sh — Downloads TCGA clinical matrices (survival, age, sex) for all 27 cancer types from UCSC Xena Hub.
- 02_download_gene_expression.sh — Downloads RNA-Seq gene expression data (HiSeqV2, log2(count+1) normalized) from UCSC Xena Hub.
- 03_download_methylation.sh — Downloads Illumina HumanMethylation450K beta values from UCSC Xena Hub.
- 04_download_mirna_expression.sh — Downloads miRNA HiSeq gene-level expression data from UCSC Xena Hub.
- 11_download_immune_signatures.sh — Downloads 68 pan-cancer immune infiltration signatures for 10,852 whitelist TCGA samples from Pan-Cancer Atlas Hub.
- 12_download_paradigm_pathway.sh — Downloads PARADIGM ssGSEA 1,387 pathway activity scores (Z-normalized) per cancer type from UCSC Xena Hub.
- 13_download_vaen_drug_data.sh — Downloads VAEN-imputed drug response data (ln IC₅₀) for TCGA samples. Supports GDSC (251 drugs, default) or CCLE (24 drugs) via -d CCLE.

### Analysis Scripts

These scripts run the core eQTM pipeline and downstream association analyses. They are located in scripts/analysis/ and can be run all at once using scripts/run_analysis.sh.

- 05_run_eqtm_pipeline.sh — Core eQTM pipeline. Runs Torch-eCpG and MatrixEQTL, applies FDR correction (BH, ≤ 0.05), filters by |r| > 0.3, annotates genomic regions (promoter / gene body / distal cis), and builds consensus CpG–gene output.
- 06_meth_drug_correlation.R — Tests CpG methylation vs. VAEN-imputed drug sensitivity (IC₅₀) using MatrixEQTL. Filter: FDR ≤ 0.05, |r| > 0.6.
- 07_immune_score_analysis.R — Tests CpG methylation vs. 68 immune infiltration signatures using MatrixEQTL. Filter: FDR ≤ 0.05, |r| > 0.6. Annotates probes with hg19 gene symbols.
- 08_mirna_probe_extraction.R — Extracts Illumina 450K methylation probes annotated to miRNA genomic regions (containing "MIR" in UCSC_RefGene_Name). Must run before script 09.
- 09_mirna_eqtl_analysis.R — Runs cis-eQTL analysis (MatrixEQTL, ±1 Mb) between miRNA-annotated CpG probes and miRNA expression. Annotates mature miRNA targets from miRBase GFF3.
- 10_paradigm_analysis.R — Tests CpG methylation vs. PARADIGM ssGSEA pathway scores using MatrixEQTL. Filter: FDR ≤ 0.05, |r| > 0.6. Annotates probes with hg19 coordinates and gene symbols.

---

## Cancer Types (27 TCGA Cohorts)

ACC, BLCA, BRCA, CHOL, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, MESO, PAAD, PCPG, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UVM

---

## Dependencies

### System Requirements

- bash >= 4.0
- curl and wget
- python3 >= 3.8
- Rscript >= 4.0
- gzip

### External Tool

Torch-eCpG — must be installed separately and placed at ./torch-ecpg/

Publication: Kober et al. (2024) BMC Bioinformatics 25:71. https://doi.org/10.1186/S12859-024-05670-4

### R Packages

Install all required packages at once:

```r
source("install_packages.R")
```

Alternatively, the same installer is also available at environment/required_packages.R.

- data.table 1.18.4 — CRAN
- MatrixEQTL 2.1.1 — CRAN
- dplyr 1.1.4 — CRAN
- RSQLite 3.53.1 — CRAN
- DBI 1.3.0 — CRAN
- matrixStats 1.5.0 — CRAN
- IlluminaHumanMethylation450kanno.ilmn12.hg19 0.6.1 — Bioconductor
- minfi 1.56.0 — Bioconductor

---

## Data Sources

- DNA methylation (450K beta values) — Illumina HumanMethylation450 — UCSC Xena Hub (TCGA)
- Gene expression (HiSeqV2, log2) — Illumina HiSeqV2 — UCSC Xena Hub (TCGA)
- miRNA expression — Illumina HiSeq — UCSC Xena Hub (TCGA)
- Clinical data (age, sex, survival) — UCSC Xena Hub (TCGA)
- Drug sensitivity (ln IC₅₀, 251 drugs) — GDSC-trained VAEN model — Jia et al. (2021) Nat Commun
- Pathway activity (1,387 pathways) — PARADIGM ssGSEA — UCSC Xena Hub (PanCan33)
- Immune infiltration (68 signatures) — Pan-Cancer Immune Landscape — Thorsson et al. (2018) Immunity
- miRNA genomic annotation — miRBase GFF3 (hsa.gff3) — https://www.mirbase.org
- Illumina 450K probe manifest — HumanMethylation450 v1.2 — Illumina / GEO

---

## Directory Structure

Set up the following layout before running. Scripts are organised into two subfolders — download/ and analysis/ — with master runner scripts at the top level.

```
project_root/
│
├── README.md
├── LICENSE
├── .zenodo.json
├── .gitignore
├── requirements.txt                        Python dependencies
├── install_packages.R                      installs all R packages in one step
│
├── scripts/
│   ├── run_downloads.sh                    runs all download scripts in order
│   ├── run_analysis.sh                     runs all analysis scripts in order
│   │
│   ├── download/                           data download scripts
│   │   ├── 01_download_clinical_data.sh
│   │   ├── 02_download_gene_expression.sh
│   │   ├── 03_download_methylation.sh
│   │   ├── 04_download_mirna_expression.sh
│   │   ├── 11_download_immune_signatures.sh
│   │   ├── 12_download_paradigm_pathway.sh
│   │   └── 13_download_vaen_drug_data.sh
│   │
│   └── analysis/                           pipeline and downstream scripts
│       ├── 05_run_eqtm_pipeline.sh
│       ├── 06_meth_drug_correlation.R
│       ├── 07_immune_score_analysis.R
│       ├── 08_mirna_probe_extraction.R
│       ├── 09_mirna_eqtl_analysis.R
│       └── 10_paradigm_analysis.R
│
├── docs/
│   └── parameter_descriptions.md          per-script parameter reference
│
├── environment/
│   └── required_packages.R                alternative R package installer
│
├── torch-ecpg/                             Torch-eCpG installation (not in zip — install separately)
│   └── Analysis/
│       └── annot/
│           ├── G.bed6                      gene genomic coordinates
│           └── M.bed6                      CpG probe genomic coordinates
│
└── data/                                   input data and annotation files (not in zip — downloaded by scripts)
    ├── ACC/
    │   ├── M.csv                           methylation matrix (probes x samples)
    │   ├── G.csv                           gene expression matrix (genes x samples)
    │   └── C.csv                           covariate matrix (age, sex, etc.)
    ├── BLCA/
    ├── BRCA/
    ├── ... (one folder per cancer type)
    ├── hsa.gff3                            miRBase annotation file (for script 09)
    ├── geneloc_miRNA.txt                   miRNA genomic coordinates (for script 09)
    ├── snpsloc.txt                         CpG probe coordinates (for script 09)
    └── humanmethylation450_15017482_v1-2.csv    Illumina 450K manifest (for script 08)
```

---

## Usage

### Step 1 — Install R Packages

```r
source("install_packages.R")
```

### Step 2 — Download All Input Data

```bash
bash scripts/run_downloads.sh
```

Or run individual scripts from scripts/download/:

```bash
bash scripts/download/01_download_clinical_data.sh
bash scripts/download/02_download_gene_expression.sh
bash scripts/download/03_download_methylation.sh
bash scripts/download/04_download_mirna_expression.sh
bash scripts/download/11_download_immune_signatures.sh
bash scripts/download/12_download_paradigm_pathway.sh
bash scripts/download/13_download_vaen_drug_data.sh          # GDSC (251 drugs, default)
# bash scripts/download/13_download_vaen_drug_data.sh -d CCLE  # CCLE (24 drugs)
```

Note: After downloading, preprocess each cancer type into M.csv (methylation), G.csv (expression), and C.csv (covariates: age, sex) inside data/{CANCER}/. Retain only primary tumor samples with complete multi-omics profiles.

### Step 3 — Run the Core cis-eQTM Pipeline

Edit DATA_DIR and TECPG_DIR in scripts/analysis/05_run_eqtm_pipeline.sh, then:

```bash
bash scripts/analysis/05_run_eqtm_pipeline.sh
```

This runs four sub-steps per cancer type:

1. Covariate cleaning — removes constant columns
2. Torch-eCpG — cis-eQTM detection within ±10 kb, FDR <= 0.05 (BH)
3. MatrixEQTL — independent validation within ±1 Mb, FDR <= 0.05 (BH)
4. Consensus — intersects both methods, filters |r| > 0.3, annotates genomic regions

### Step 4 — Run Downstream Association Analyses

```bash
bash scripts/run_analysis.sh
```

Or run individually:

```bash
Rscript scripts/analysis/06_meth_drug_correlation.R      # CpG vs drug sensitivity
Rscript scripts/analysis/07_immune_score_analysis.R      # CpG vs immune infiltration
Rscript scripts/analysis/08_mirna_probe_extraction.R     # extract miRNA-annotated probes
Rscript scripts/analysis/09_mirna_eqtl_analysis.R        # CpG vs miRNA expression (cis-eQTL)
Rscript scripts/analysis/10_paradigm_analysis.R          # CpG vs PARADIGM pathway scores
```

---

## Output Files

- 05_run_eqtm_pipeline.sh — data/{CANCER}/Results/final_{CANCER}.csv
- 06_meth_drug_correlation.R — results/drug_associations/{CANCER}_Meth_Drug_Associations.csv and .sqlite
- 07_immune_score_analysis.R — results/immune_associations/{CANCER}_Immune_Associations.csv and .sqlite
- 08_mirna_probe_extraction.R — results/mirna_probes/M_{CANCER}_miRNA_annotated.csv
- 09_mirna_eqtl_analysis.R — results/mirna_eqtl/{CANCER}/{CANCER}_cis_mirna_FDR_0.05_annotated.csv
- 10_paradigm_analysis.R — results/paradigm_associations/{CANCER}_Paradigm_MatrixEQTL_FDR0.05_Corr0.6.csv and .sqlite

### Key Columns in the Final eQTM File (final_{CANCER}.csv)

- mt_id — CpG probe ID (Illumina 450K)
- gt_id — Associated gene ID
- mt_chrom, mt_chromStart, mt_strand — Probe genomic coordinates
- gt_chrom, gt_chromStart, gt_strand — Gene genomic coordinates
- region — Genomic context: CIS / PROMOTER / DISTAL
- mt_est, mt_err, mt_t, mt_p — Torch-eCpG linear regression statistics
- beta, t-stat, p-value, FDR — MatrixEQTL statistics
- correlation — Pearson correlation coefficient (r)

---

## Statistical Filtering Criteria

All analyses apply a false discovery rate (FDR) threshold of <= 0.05 using the Benjamini–Hochberg (BH) method. For the consensus eQTM module, an additional Pearson correlation filter of |r| > 0.3 is applied. Drug sensitivity, immune infiltration, and PARADIGM pathway associations each require a stricter correlation threshold of |r| > 0.6. The miRNA cis-eQTL module (script 09) applies FDR <= 0.05 only — no additional correlation filter is used.

---

## Funding

- Department of Biotechnology (DBT), Government of India — BT/PR52559/MED/30/2534/2024
- Indian Council of Medical Research (ICMR), Government of India — IIRPSG-2024-01-02098
- B.T.K. is supported by a University Grants Commission (UGC) fellowship, Government of India

---

## License

This code is released under the MIT License — see LICENSE for details.

---

## Contact

Rahul Kumar | rahulk@bt.iith.ac.in
Department of Biotechnology, Indian Institute of Technology Hyderabad, India
