[Uploading README.md…]()
# Onco-eQTM: Pan-Cancer Expression Quantitative Trait Methylation Pipeline


---

## Overview

This repository contains all analysis scripts used to build the **Onco-eQTM** database — a pan-cancer resource that maps cis-regulatory DNA methylation to gene expression, drug sensitivity, immune infiltration, miRNA regulation, and pathway activity across **27 TCGA cancer types** and **6,880 primary tumor samples**.

The pipeline identifies high-confidence CpG–gene associations (eQTMs) using a consensus strategy combining **Torch-eCpG** and **MatrixEQTL**, then integrates results with four additional functional layers:

| Layer | Scale |
|-------|-------|
| cis-eQTM associations (CpG ↔ gene expression) | 5.25 million |
| Drug response associations (251 compounds, GDSC) | 14.45 million |
| Pathway activity associations (1,387 PARADIGM pathways) | 13.6 million |
| Immune infiltration associations (68 immune cell types) | 3.55 million |
| miRNA-associated CpG associations | 4.52 million |

The web database is freely available at: **https://project.iith.ac.in/cgntlab/OncoeQTM/**

---

## Citation

If you use these scripts or data, please cite:

> Korra BT, Nishana M, Kumar R. *Expression Quantitative Trait Methylation Across Multiple Cancer Types with Functional and Therapeutic Characterization using Onco-eQTM*. [Journal, Year, DOI — to be updated upon publication]

---

## Repository Contents

### Download Scripts (Bash)

| Script | Description |
|--------|-------------|
| `01_download_clinical_data.sh` | Downloads TCGA clinical matrices (survival, age, sex) for all 27 cancer types from UCSC Xena Hub |
| `02_download_gene_expression.sh` | Downloads RNA-Seq gene expression data (HiSeqV2, log2(count+1) normalized) from UCSC Xena Hub |
| `03_download_methylation.sh` | Downloads Illumina HumanMethylation450K beta values from UCSC Xena Hub |
| `04_download_mirna_expression.sh` | Downloads miRNA HiSeq gene-level expression data from UCSC Xena Hub |
| `11_download_immune_signatures.sh` | Downloads 68 pan-cancer immune infiltration signatures for 10,852 whitelist TCGA samples from Pan-Cancer Atlas Hub |
| `12_download_paradigm_pathway.sh` | Downloads PARADIGM ssGSEA 1,387 pathway activity scores (Z-normalized) per cancer type from UCSC Xena Hub |
| `13_download_vaen_drug_data.sh` | Downloads VAEN-imputed drug response data (ln IC₅₀) for TCGA samples; supports GDSC (251 drugs, default) or CCLE (24 drugs) |

### Analysis Scripts

| Script | Description |
|--------|-------------|
| `05_run_eqtm_pipeline.sh` | Core eQTM pipeline: runs Torch-eCpG and MatrixEQTL, applies FDR correction (BH, ≤0.05), filters by \|r\| > 0.3, annotates genomic regions (promoter / gene body / distal cis), and builds consensus CpG–gene output |
| `06_meth_drug_correlation.R` | Tests CpG methylation vs. VAEN-imputed drug sensitivity (IC₅₀) using MatrixEQTL; filters FDR ≤ 0.05, \|r\| > 0.6 |
| `07_immune_score_analysis.R` | Tests CpG methylation vs. 68 immune infiltration signatures using MatrixEQTL; filters FDR ≤ 0.05, \|r\| > 0.6; annotates probes with hg19 gene symbols |
| `08_mirna_probe_extraction.R` | Extracts Illumina 450K methylation probes annotated to miRNA genomic regions (containing "MIR" in UCSC_RefGene_Name) |
| `09_mirna_eqtl_analysis.R` | Runs cis-eQTL analysis (MatrixEQTL, ±1 Mb) between miRNA-annotated CpG probes and miRNA expression; annotates mature miRNA targets from miRBase GFF3 |
| `10_paradigm_analysis.R` | Tests CpG methylation vs. PARADIGM ssGSEA pathway scores using MatrixEQTL; filters FDR ≤ 0.05, \|r\| > 0.6; annotates probes with hg19 coordinates and gene symbols |
| `required_packages.R` | Installs all required R packages (CRAN + Bioconductor) in a single step |

---

## Cancer Types (27 TCGA Cohorts)

ACC, BLCA, BRCA, CHOL, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, MESO, PAAD, PCPG, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UVM

---

## Dependencies

### System Requirements
- `bash` ≥ 4.0
- `curl` and `wget`
- `python3` ≥ 3.8
- `Rscript` ≥ 4.0
- `gzip`

### External Tool
- **Torch-eCpG** — must be installed separately and placed at `./torch-ecpg/`
  Publication: Kober et al. (2024) *BMC Bioinformatics* 25:71. https://doi.org/10.1186/S12859-024-05670-4

### R Packages
Install all required packages at once:

```r
source("required_packages.R")
```

Packages installed:

| Package | Version | Source |
|---------|---------|--------|
| `data.table` | 1.18.4 | CRAN |
| `MatrixEQTL` | 2.1.1 | CRAN |
| `dplyr` | 1.1.4 | CRAN |
| `RSQLite` | 3.53.1 | CRAN |
| `DBI` | 1.3.0 | CRAN |
| `matrixStats` | 1.5.0 | CRAN |
| `IlluminaHumanMethylation450kanno.ilmn12.hg19` | 0.6.1 | Bioconductor |
| `minfi` | 1.56.0 | Bioconductor |

---

## Data Sources

| Data | Platform | Source |
|------|----------|--------|
| DNA methylation (450K beta values) | Illumina HumanMethylation450 | UCSC Xena Hub (TCGA) |
| Gene expression (HiSeqV2, log2) | Illumina HiSeqV2 | UCSC Xena Hub (TCGA) |
| miRNA expression | Illumina HiSeq | UCSC Xena Hub (TCGA) |
| Clinical data (age, sex, survival) | — | UCSC Xena Hub (TCGA) |
| Drug sensitivity (ln IC₅₀, 251 drugs) | GDSC-trained VAEN model | Jia et al. (2021) *Nat Commun* |
| Pathway activity (1,387 pathways) | PARADIGM ssGSEA | UCSC Xena Hub (PanCan33) |
| Immune infiltration (68 signatures) | Pan-Cancer Immune Landscape | Thorsson et al. (2018) *Immunity* |
| miRNA genomic annotation | miRBase GFF3 (hsa.gff3) | https://www.mirbase.org |
| Illumina 450K probe manifest | HumanMethylation450 v1.2 | Illumina / GEO |

---

## Directory Structure

Set up the following layout before running:

```
project_root/
├── scripts/                          ← this repository
│   ├── 01_download_clinical_data.sh
│   ├── 02_download_gene_expression.sh
│   ├── 03_download_methylation.sh
│   ├── 04_download_mirna_expression.sh
│   ├── 05_run_eqtm_pipeline.sh
│   ├── 06_meth_drug_correlation.R
│   ├── 07_immune_score_analysis.R
│   ├── 08_mirna_probe_extraction.R
│   ├── 09_mirna_eqtl_analysis.R
│   ├── 10_paradigm_analysis.R
│   ├── 11_download_immune_signatures.sh
│   ├── 12_download_paradigm_pathway.sh
│   ├── 13_download_vaen_drug_data.sh
│   └── required_packages.R
├── torch-ecpg/                       ← Torch-eCpG installation
│   └── Analysis/
│       └── annot/
│           ├── G.bed6                ← gene genomic coordinates
│           └── M.bed6                ← CpG probe genomic coordinates
└── data/
    ├── ACC/
    │   ├── M.csv                     ← methylation matrix  (probes × samples)
    │   ├── G.csv                     ← gene expression matrix (genes × samples)
    │   └── C.csv                     ← covariate matrix (samples: age, sex)
    ├── BLCA/ ...
    ├── [one folder per cancer type]
    ├── hsa.gff3                      ← miRBase annotation (for script 09)
    ├── geneloc_miRNA.txt             ← miRNA genomic coordinates (for script 09)
    ├── snpsloc.txt                   ← CpG probe coordinates (for script 09)
    └── humanmethylation450_15017482_v1-2.csv   ← Illumina manifest (for script 08)
```

---

## Usage

### Step 1 — Install R packages

```r
source("required_packages.R")
```

### Step 2 — Download all raw data

```bash
bash 01_download_clinical_data.sh
bash 02_download_gene_expression.sh
bash 03_download_methylation.sh
bash 04_download_mirna_expression.sh
bash 11_download_immune_signatures.sh
bash 12_download_paradigm_pathway.sh
bash 13_download_vaen_drug_data.sh          # GDSC (251 drugs, default)
# bash 13_download_vaen_drug_data.sh -d CCLE  # CCLE (24 drugs)
```

> **Note:** After downloading, preprocess each cancer type's data into `M.csv` (methylation), `G.csv` (expression), and `C.csv` (covariates: age, sex) inside `data/{CANCER}/`. Only primary tumor samples with complete multi-omics profiles are retained. Clinical covariates are limited to age and sex.

### Step 3 — Run the core eQTM pipeline

Edit `DATA_DIR` and `TECPG_DIR` in `05_run_eqtm_pipeline.sh`, then:

```bash
bash 05_run_eqtm_pipeline.sh
```

This runs four sub-steps per cancer type:
1. **Covariate cleaning** — removes constant columns
2. **Torch-eCpG** — cis-eQTM detection within ±10 kb; FDR ≤ 0.05 (BH)
3. **MatrixEQTL** — independent validation within ±1 Mb; FDR ≤ 0.05 (BH)
4. **Consensus** — intersects both methods; filters \|r\| > 0.3; annotates genomic regions

### Step 4 — Run downstream association analyses

```bash
Rscript 06_meth_drug_correlation.R      # CpG ↔ drug sensitivity
Rscript 07_immune_score_analysis.R      # CpG ↔ immune infiltration
Rscript 08_mirna_probe_extraction.R     # extract miRNA-annotated probes
Rscript 09_mirna_eqtl_analysis.R        # CpG ↔ miRNA expression (cis-eQTL)
Rscript 10_paradigm_analysis.R          # CpG ↔ PARADIGM pathway scores
```

---

## Output Files

| Script | Output Path | Format |
|--------|-------------|--------|
| `05_run_eqtm_pipeline.sh` | `data/{CANCER}/Results/final_{CANCER}.csv` | CSV |
| `06_meth_drug_correlation.R` | `results/drug_associations/{CANCER}_Meth_Drug_Associations.csv/.sqlite` | CSV + SQLite |
| `07_immune_score_analysis.R` | `results/immune_associations/{CANCER}_Immune_Associations.csv/.sqlite` | CSV + SQLite |
| `08_mirna_probe_extraction.R` | `results/mirna_probes/M_{CANCER}_miRNA_annotated.csv` | CSV |
| `09_mirna_eqtl_analysis.R` | `results/mirna_eqtl/{CANCER}/{CANCER}_cis_mirna_FDR_0.05_annotated.csv` | CSV |
| `10_paradigm_analysis.R` | `results/paradigm_associations/{CANCER}_Paradigm_MatrixEQTL_FDR0.05_Corr0.6.csv/.sqlite` | CSV + SQLite |

### Key columns in the final eQTM file (`final_{CANCER}.csv`)

| Column | Description |
|--------|-------------|
| `mt_id` | CpG probe ID (Illumina 450K) |
| `gt_id` | Associated gene ID |
| `mt_chrom`, `mt_chromStart`, `mt_strand` | Probe genomic coordinates |
| `gt_chrom`, `gt_chromStart`, `gt_strand` | Gene genomic coordinates |
| `region` | Genomic context: CIS / PROMOTER / DISTAL |
| `mt_est`, `mt_err`, `mt_t`, `mt_p` | Torch-eCpG linear regression statistics |
| `beta`, `t-stat`, `p-value`, `FDR` | MatrixEQTL statistics |
| `correlation` | Pearson correlation coefficient (r) |

---

## Statistical Filtering Criteria

All analyses apply a false discovery rate (FDR) threshold of ≤ 0.05 using the Benjamini–Hochberg (BH) multiple-testing correction method throughout the workflow. For the consensus eQTM module, an additional Pearson correlation filter of |r| > 0.3 is applied. Drug sensitivity, immune infiltration, PARADIGM pathway, and miRNA cis-eQTL associations each require a stricter correlation threshold of |r| > 0.6.

---

## Funding

- Department of Biotechnology (DBT), Government of India — BT/PR52559/MED/30/2534/2024
- Indian Council of Medical Research (ICMR), Government of India — IIRPSG-2024-01-02098
- B.T.K. is supported by a University Grants Commission (UGC) fellowship, Government of India

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

**Rahul Kumar** | rahulk@bt.iith.ac.in
Department of Biotechnology, Indian Institute of Technology Hyderabad, India
