# Parameter Descriptions

## Core eQTM Pipeline (`05_run_eqtm_pipeline.sh`)

| Parameter          | Default       | Description                                 |
| ------------------ | ------------- | ------------------------------------------- |
| `DATA_DIR`         | `./data`      | Root directory containing cancer subfolders |
| `TECPG_DIR`        | `./torch-ecpg`| Path to Torch-eCpG installation             |
| `cisDist`          | `1,000,000`   | Cis window in base pairs (1 Mb)             |
| `cpu-threads`      | `16`          | Number of CPU threads for Torch-eCpG        |
| FDR cutoff         | `0.05`        | Benjamini-Hochberg threshold                |
| Correlation cutoff | `0.3`         | Minimum absolute Pearson correlation        |

## Drug Sensitivity (`06_meth_drug_correlation.R`)

| Parameter     | Default                       | Description                              |
| ------------- | ----------------------------- | ---------------------------------------- |
| `FDR_CUTOFF`  | `0.05`                        | BH-corrected p-value threshold           |
| `CORR_CUTOFF` | `0.6`                         | Minimum absolute correlation             |
| `DRUG_DIR`    | `./data`                      | Folder with `CANCER_drug_data.csv` files |
| `OUTPUT_DIR`  | `./results/drug_associations` | Output location                          |

## Immune Infiltration (`07_immune_score_analysis.R`)

| Parameter     | Default                         | Description                                     |
| ------------- | ------------------------------- | ----------------------------------------------- |
| `FDR_CUTOFF`  | `0.05`                          | BH-corrected p-value threshold                  |
| `CORR_CUTOFF` | `0.6`                           | Minimum absolute correlation                      |
| `OUTPUT_DIR`  | `./results/immune_associations` | Output location                                 |

## PARADIGM Pathways (`10_paradigm_analysis.R`)

| Parameter     | Default                           | Description                    |
| ------------- | --------------------------------- | ------------------------------ |
| `FDR_CUTOFF`  | `0.05`                            | BH-corrected p-value threshold |
| `CORR_CUTOFF` | `0.6`                             | Minimum absolute correlation   |
| `OUTPUT_DIR`  | `./results/paradigm_associations` | Output location                |
