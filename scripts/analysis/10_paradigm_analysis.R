#!/usr/bin/env Rscript
# =============================================================================
# Pan-Cancer Methylation vs PARADIGM Pathway Analysis
#
# Loops through all 27 TCGA cancer types to run MatrixEQTL testing
# the correlation between significant eQTM CpG probes and PARADIGM pathway
# activities (ssGSEA Z-scores).
#
# Filters results for FDR < 0.05 and |correlation| > 0.6.
# Annotates probes with genomic coordinates and gene symbols using hg19.
#
# Output: CSV + SQLite database in results/paradigm_associations/
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

# Base directory where your TCGA cancer folders are
DATA_DIR     <- "./data"

# Directory where output results will be stored
MAIN_OUT_DIR <- "./results/paradigm_associations"

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

if (!dir.exists(MAIN_OUT_DIR)) dir.create(MAIN_OUT_DIR, recursive = TRUE)

# --- Load hg19 Probe Annotation (runs once) ---
cat("Loading hg19 annotation data...\n")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_dt  <- data.table(
  mt_id       = rownames(ann450k),
  chr         = ann450k$chr,
  pos         = ann450k$pos,
  gene_symbol = ann450k$UCSC_RefGene_Name
)
ann_dt$gene_symbol <- sapply(strsplit(ann_dt$gene_symbol, ";"), function(x) x[1])
ann_dt$gene_symbol[ann_dt$gene_symbol == ""] <- NA

cat("\n=== Pan-Cancer Methylation-PARADIGM Analysis ===\n\n")

for (cancer_type in CANCERS) {
  
  tryCatch({
    cat("Processing:", cancer_type, "\n")
    
    # --- Define Paths ---
    cancer_dir   <- file.path(DATA_DIR, cancer_type)
    input_file   <- file.path(cancer_dir, "Results", paste0("final_", cancer_type, ".csv"))
    meth_file    <- file.path(cancer_dir, "M.csv")
    cov_file     <- file.path(cancer_dir, "C.csv")
    
    # ssGSEA Paradigm pathway activity file matching the cancer code
    paradigm_filename <- paste0(cancer_type, "_PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level_Z.txt.gz")
    paradigm_file     <- file.path(cancer_dir, paradigm_filename)
    
    # Validate Files
    if (!all(file.exists(input_file, meth_file, cov_file, paradigm_file))) {
      cat("  Skipping:", cancer_type, "- required files missing.\n")
      next
    }
    
    # --- Load Data ---
    merged_data      <- fread(input_file)
    unique_probes    <- unique(merged_data$mt_id)
    
    methylation_data <- fread(meth_file)
    if (colnames(methylation_data)[1] != "sample") {
      colnames(methylation_data)[1] <- "sample"
    }
    filtered_meth <- methylation_data[sample %in% unique_probes, ]
    rm(methylation_data); gc()
    
    # Clean methylation sample IDs
    meth_samples <- colnames(filtered_meth)[-1]
    meth_samples_clean <- gsub("\\.", "-", meth_samples)
    colnames(filtered_meth) <- c("sample", meth_samples_clean)
    
    # Load and clean paradigm pathway file
    paradigm_raw <- fread(paradigm_file)
    paradigm_df  <- as.data.frame(paradigm_raw)
    rownames(paradigm_df) <- paradigm_df[, 1]
    paradigm_df  <- paradigm_df[, -1]
    
    para_samples_clean <- gsub("\\.", "-", colnames(paradigm_df))
    colnames(paradigm_df) <- para_samples_clean
    
    # Load covariates
    covariates <- fread(cov_file)
    colnames(covariates)[1] <- "sampleID"
    covariates$sampleID <- gsub("\\.", "-", covariates$sampleID)
    
    # --- Align Samples ---
    common_samples <- Reduce(intersect, list(
      meth_samples_clean,
      para_samples_clean,
      covariates$sampleID
    ))
    
    if (length(common_samples) < 10) {
      cat("  Too few common samples ( < 10 ) - skipping.\n")
      next
    }
    
    cat("  Probes:", length(unique_probes), " | Samples:", length(common_samples), "\n")
    
    # Subset and format matrices
    meth_df <- as.data.frame(filtered_meth)
    rownames(meth_df) <- meth_df$sample; meth_df$sample <- NULL
    meth_matrix <- as.matrix(meth_df)[, common_samples, drop = FALSE]
    
    paradigm_matrix <- as.matrix(paradigm_df)[, common_samples, drop = FALSE]
    
    cov_sub <- covariates[match(common_samples, covariates$sampleID)]
    cov_matrix <- as.matrix(t(cov_sub[, -"sampleID", with = FALSE]))
    colnames(cov_matrix) <- common_samples
    
    # --- Run MatrixEQTL ---
    cat("  Running MatrixEQTL...\n")
    snpspos <- data.frame(SNP = rownames(meth_matrix), chr = 1, pos = 1:nrow(meth_matrix))
    genepos <- data.frame(gene = rownames(paradigm_matrix), chr = 1, pos = 1:nrow(paradigm_matrix))
    
    me <- Matrix_eQTL_main(
      snps = SlicedData$new(meth_matrix),
      gene = SlicedData$new(paradigm_matrix),
      cvrt = SlicedData$new(cov_matrix),
      output_file_name = "",
      pvOutputThreshold = 1,
      useModel = modelLINEAR,
      verbose = FALSE, pvalue.hist = FALSE
    )
    
    # --- Process and Filter Results ---
    results <- me$all$eqtls
    N <- length(common_samples)
    K <- nrow(cov_matrix)
    results$correlation <- results$statistic / sqrt(results$statistic^2 + N - K - 2)
    results$SE <- results$beta / results$statistic
    
    sig_results <- results[results$FDR < FDR_CUTOFF & abs(results$correlation) > CORR_CUTOFF, ]
    sig_results <- sig_results[order(sig_results$FDR), ]
    
    # --- Annotate and Save ---
    if (nrow(sig_results) > 0) {
      sig_anno <- sig_results %>%
        left_join(ann_dt, by = c("snps" = "mt_id")) %>%
        rename(
          CpG_ID = snps, Pathway = gene, 
          Methylation_Chr = chr, Methylation_Pos = pos, Associated_Gene = gene_symbol
        ) %>%
        select(CpG_ID, Methylation_Chr, Methylation_Pos, Associated_Gene, Pathway, 
               beta, statistic, pvalue, FDR, correlation, SE)
      sig_anno$Cancer_Type <- cancer_type
      
      out_base <- paste0(cancer_type, "_Paradigm_MatrixEQTL_FDR", FDR_CUTOFF, "_Corr", CORR_CUTOFF)
      out_csv  <- file.path(MAIN_OUT_DIR, paste0(out_base, ".csv"))
      out_db   <- file.path(MAIN_OUT_DIR, paste0(out_base, ".sqlite"))
      
      fwrite(sig_anno, out_csv)
      
      con <- dbConnect(RSQLite::SQLite(), out_db)
      dbWriteTable(con, "paradigm_results", sig_anno, overwrite = TRUE)
      dbDisconnect(con)
      
      cat("  Saved:", nrow(sig_anno), "associations.\n")
    } else {
      cat("  No significant associations found.\n")
    }
    
  }, error = function(e) {
    cat("  ERROR in", cancer_type, ":", conditionMessage(e), "\n")
  })
  
  gc()
}

cat("\n=== PARADIGM Analysis Complete ===\n")
