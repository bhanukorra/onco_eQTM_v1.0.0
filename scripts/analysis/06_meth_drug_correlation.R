#!/usr/bin/env Rscript
# =============================================================================
# Pan-Cancer Methylation vs Drug Sensitivity Analysis
#
# For each cancer, this script tests whether the methylation level of
# significant eQTM CpG probes is correlated with drug sensitivity (IC50).
#
# It uses MatrixEQTL to run the association test, then filters for
# FDR < 0.05 and |correlation| > 0.6.
#
# Output: CSV + SQLite per cancer type.
# =============================================================================

library(data.table)
library(MatrixEQTL)
library(dplyr)
library(RSQLite)
library(DBI)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(normalizePath(file.path(script_dir, "..", "..")))
}

# =============================================================================
# SETTINGS — change these paths to match your system
# =============================================================================

# Folder containing your cancer data (each subfolder has M.csv, C.csv, Results/)
DATA_DIR     <- "./data"

# Folder containing drug sensitivity files (CANCER_drug_data.csv per cancer)
DRUG_DIR     <- "./data"

# Where to save the drug association results
OUTPUT_DIR   <- "./results/drug_associations"

# =============================================================================
# Parameters
# =============================================================================
FDR_CUTOFF   <- 0.05
CORR_CUTOFF  <- 0.6

# 27 TCGA cancer types
CANCERS <- c(
  "ACC", "BLCA", "BRCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UVM"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("\n=== Pan-Cancer Methylation-Drug Analysis ===\n\n")

for (cancer in CANCERS) {

  cat("Processing:", cancer, "\n")

  # --- Locate files ---
  results_dir <- file.path(DATA_DIR, cancer, "Results")
  final_file  <- file.path(results_dir, paste0("final_", cancer, ".csv"))
  meth_file   <- file.path(DATA_DIR, cancer, "M.csv")
  cov_file    <- file.path(DATA_DIR, cancer, "C.csv")
  drug_file   <- file.path(DRUG_DIR, cancer, paste0(cancer, "_drug_data.csv"))

  # Skip if any file is missing
  if (!all(file.exists(final_file, meth_file, cov_file, drug_file))) {
    cat("  Skipping:", cancer, "- required files missing.\n")
    next
  }

  # --- Load significant eQTM probes ---
  final_res    <- fread(final_file, select = c("mt_id", "FDR"))
  unique_probes <- unique(final_res[FDR < 0.05]$mt_id)
  rm(final_res); gc()

  if (length(unique_probes) == 0) {
    cat("  No significant eQTM probes for", cancer, "- skipping.\n")
    next
  }

  # --- Load and filter methylation data ---
  meth_data <- fread(meth_file)
  setnames(meth_data, colnames(meth_data)[1], "probe")
  meth_filtered <- meth_data[probe %in% unique_probes]
  rm(meth_data); gc()

  # --- Load drug and covariate data ---
  drug_data <- fread(drug_file)
  setnames(drug_data, colnames(drug_data)[1], "sampleID")

  covariates <- fread(cov_file)
  if (!("sampleID" %in% colnames(covariates))) {
    setnames(covariates, colnames(covariates)[1], "sampleID")
  }

  # --- Find common samples across all three datasets ---
  common_samples <- Reduce(intersect, list(
    colnames(meth_filtered)[-1],
    drug_data$sampleID,
    covariates$sampleID
  ))

  if (length(common_samples) < 10) {
    cat("  Too few common samples (", length(common_samples), ") - skipping.\n")
    next
  }

  cat("  Probes:", length(unique_probes), " | Samples:", length(common_samples), "\n")

  # --- Build matrices for MatrixEQTL ---
  meth_mat <- as.matrix(meth_filtered[, ..common_samples])
  rownames(meth_mat) <- meth_filtered$probe

  drug_sub <- drug_data[match(common_samples, drug_data$sampleID)]
  drug_mat <- as.matrix(t(drug_sub[, -"sampleID", with = FALSE]))
  colnames(drug_mat) <- common_samples

  cov_sub  <- covariates[match(common_samples, covariates$sampleID)]
  cov_mat  <- as.matrix(t(cov_sub[, -"sampleID", with = FALSE]))
  colnames(cov_mat) <- common_samples

  # --- Run MatrixEQTL ---
  cat("  Running MatrixEQTL...\n")
  me <- Matrix_eQTL_main(
    snps = SlicedData$new(meth_mat),
    gene = SlicedData$new(drug_mat),
    cvrt = SlicedData$new(cov_mat),
    output_file_name = NULL,
    pvOutputThreshold = 0.05,
    useModel = modelLINEAR,
    verbose = FALSE,
    pvalue.hist = FALSE
  )

  results <- me$all$eqtls
  if (nrow(results) == 0) { cat("  No associations found.\n"); next }

  # --- Compute correlation and filter ---
  df_res <- length(common_samples) - nrow(cov_mat) - 2
  results$correlation <- results$statistic / sqrt(results$statistic^2 + df_res)
  results_filtered <- results[results$FDR < FDR_CUTOFF & abs(results$correlation) > CORR_CUTOFF, ]

  if (nrow(results_filtered) == 0) { cat("  No associations passed filters.\n"); next }

  # --- Save CSV ---
  out_csv <- file.path(OUTPUT_DIR, paste0(cancer, "_Meth_Drug_Associations.csv"))
  fwrite(results_filtered, out_csv)

  # --- Save SQLite ---
  out_db <- file.path(OUTPUT_DIR, paste0(cancer, "_Meth_Drug_Associations.sqlite"))
  con <- dbConnect(RSQLite::SQLite(), out_db)
  dbWriteTable(con, "meth_drug", results_filtered, overwrite = TRUE)
  dbDisconnect(con)

  cat("  Saved:", nrow(results_filtered), "associations.\n")
  rm(meth_filtered, drug_data, covariates, meth_mat, drug_mat, cov_mat, results, results_filtered); gc()
}

cat("\n=== Drug Analysis Complete ===\n")
