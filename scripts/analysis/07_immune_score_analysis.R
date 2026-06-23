#!/usr/bin/env Rscript
# =============================================================================
# Pan-Cancer Methylation vs Immune Infiltration Analysis
#
# For each cancer, this script tests whether the methylation level of
# significant eQTM CpG probes is correlated with immune cell infiltration
# scores (68 immune signatures from the Pan-Cancer Immune Landscape).
#
# It uses MatrixEQTL for association testing, filters for FDR < 0.05 and
# |correlation| > 0.6, and annotates results with hg19 gene names.
#
# Output: CSV + SQLite per cancer type.
# =============================================================================

library(data.table)
library(MatrixEQTL)
library(matrixStats)
library(RSQLite)
library(DBI)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(normalizePath(file.path(script_dir, "..", "..")))
}

# =============================================================================
# SETTINGS — change these paths to match your system
# =============================================================================

# Folder containing cancer data (each subfolder has M.csv, C.csv, etc.)
DATA_DIR   <- "./data"

# Where to save immune association results
OUTPUT_DIR <- "./results/immune_associations"

# =============================================================================
# Parameters
# =============================================================================
FDR_CUTOFF  <- 0.05
CORR_CUTOFF <- 0.6

# 27 TCGA cancer types
CANCERS <- c(
  "ACC", "BLCA", "BRCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UVM"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- Load hg19 Probe Annotation (runs once) ---
cat("Loading hg19 annotation...\n")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_dt  <- data.table(
  mt_id       = rownames(ann450k),
  chr         = ann450k$chr,
  pos         = ann450k$pos,
  gene_symbol = ann450k$UCSC_RefGene_Name
)
ann_dt$gene_symbol <- sapply(strsplit(ann_dt$gene_symbol, ";"), function(x) x[1])
ann_dt$gene_symbol[ann_dt$gene_symbol == ""] <- NA

cat("\n=== Pan-Cancer Methylation-Immune Analysis ===\n\n")

for (cancer in CANCERS) {

  tryCatch({
    cat("Processing:", cancer, "\n")

    # --- Locate files ---
    cancer_dir  <- file.path(DATA_DIR, cancer)
    input_file  <- file.path(cancer_dir, "Results", paste0("final_", cancer, ".csv"))
    meth_file   <- file.path(cancer_dir, "M.csv")
    cov_file    <- file.path(cancer_dir, "C.csv")
    immune_file <- file.path(cancer_dir, paste0(cancer, "_immune_signatures.csv"))

    # Skip if files missing
    if (!all(file.exists(input_file, meth_file, cov_file, immune_file))) {
      cat("  Skipping:", cancer, "- files missing.\n"); next
    }

    # --- Load eQTM probes ---
    merged_data   <- fread(input_file)
    unique_probes <- unique(merged_data$mt_id)

    # --- Load methylation data ---
    meth_data <- fread(meth_file)
    if (colnames(meth_data)[1] != "sample") colnames(meth_data)[1] <- "sample"
    meth_filtered <- meth_data[sample %in% unique_probes, ]
    rm(meth_data); gc()

    # Clean sample IDs (dots to dashes)
    meth_samples <- colnames(meth_filtered)[-1]
    colnames(meth_filtered) <- c("sample", gsub("\\.", "-", meth_samples))

    # --- Load and transpose immune data ---
    immune_raw <- fread(immune_file)
    immune_df  <- as.data.frame(immune_raw)
    rownames(immune_df) <- immune_df[, 1]
    immune_df  <- immune_df[, -1]
    immune_matrix <- t(immune_df)
    colnames(immune_matrix) <- gsub("\\.", "-", colnames(immune_matrix))

    # --- Load covariates ---
    covariates <- fread(cov_file)
    colnames(covariates)[1] <- "sampleID"
    covariates$sampleID <- gsub("\\.", "-", covariates$sampleID)

    # --- Find common samples ---
    common <- Reduce(intersect, list(
      colnames(meth_filtered)[-1],
      colnames(immune_matrix),
      covariates$sampleID
    ))

    if (length(common) < 10) { cat("  Too few samples - skipping.\n"); next }
    cat("  Probes:", length(unique_probes), " | Samples:", length(common), "\n")

    # --- Build matrices ---
    meth_df <- as.data.frame(meth_filtered)
    rownames(meth_df) <- meth_df$sample; meth_df$sample <- NULL
    meth_mat <- as.matrix(meth_df)[, common, drop = FALSE]

    imm_mat <- immune_matrix[, common, drop = FALSE]

    cov_sub <- covariates[match(common, covariates$sampleID)]
    cov_mat <- as.matrix(t(cov_sub[, -"sampleID", with = FALSE]))
    colnames(cov_mat) <- common

    # --- Run MatrixEQTL ---
    cat("  Running MatrixEQTL...\n")
    me <- Matrix_eQTL_main(
      snps = SlicedData$new(meth_mat),
      gene = SlicedData$new(imm_mat),
      cvrt = SlicedData$new(cov_mat),
      output_file_name = "",
      pvOutputThreshold = 1,
      useModel = modelLINEAR,
      verbose = FALSE, pvalue.hist = FALSE
    )

    results <- me$all$eqtls
    N <- length(common); K <- nrow(cov_mat)
    results$correlation <- results$statistic / sqrt(results$statistic^2 + N - K - 2)

    # --- Filter and annotate ---
    sig <- results[results$FDR < FDR_CUTOFF & abs(results$correlation) > CORR_CUTOFF, ]
    sig <- sig[order(sig$FDR), ]

    if (nrow(sig) > 0) {
      sig_anno <- sig %>%
        left_join(ann_dt, by = c("snps" = "mt_id")) %>%
        rename(CpG_ID = snps, Immune_Signature = gene,
               Methylation_Chr = chr, Methylation_Pos = pos, Associated_Gene = gene_symbol) %>%
        select(CpG_ID, Methylation_Chr, Methylation_Pos, Associated_Gene,
               Immune_Signature, beta, statistic, pvalue, FDR, correlation)
      sig_anno$Cancer_Type <- cancer

      # Save CSV
      out_csv <- file.path(OUTPUT_DIR, paste0(cancer, "_Immune_Associations.csv"))
      fwrite(sig_anno, out_csv)

      # Save SQLite
      out_db <- file.path(OUTPUT_DIR, paste0(cancer, "_Immune_Associations.sqlite"))
      con <- dbConnect(RSQLite::SQLite(), out_db)
      dbWriteTable(con, "immune_results", sig_anno, overwrite = TRUE)
      dbDisconnect(con)

      cat("  Saved:", nrow(sig_anno), "associations.\n")
    } else {
      cat("  No significant associations.\n")
    }

  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })

  gc()
}

cat("\n=== Immune Analysis Complete ===\n")
