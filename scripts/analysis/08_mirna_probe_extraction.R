#!/usr/bin/env Rscript
# =============================================================================
# Pan-Cancer miRNA Probe Extraction
#
# This script extracts methylation probes (CpG sites) that are annotated to
# miRNA genes (containing "MIR" in the Illumina 450K manifest).
#
# For each cancer, it reads M.csv, filters for miRNA-annotated probes,
# and saves the filtered methylation matrix.
#
# Output: M_CANCER_miRNA_annotated.csv per cancer type.
# =============================================================================

library(data.table)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  setwd(normalizePath(file.path(script_dir, "..", "..")))
}

# =============================================================================
# SETTINGS — change these paths to match your system
# =============================================================================

# Path to Illumina HumanMethylation450 manifest file (v1.2)
MANIFEST_FILE <- "./data/humanmethylation450_15017482_v1-2.csv"

# Folder containing cancer data (each subfolder has M.csv)
DATA_DIR      <- "./data"

# Where to save the miRNA-filtered methylation files
OUTPUT_DIR    <- "./results/mirna_probes"

# =============================================================================

# 27 TCGA cancer types
CANCERS <- c(
  "ACC", "BLCA", "BRCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UVM"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- Step 1: Load manifest and identify miRNA probes ---
cat("Loading manifest:", MANIFEST_FILE, "\n")
manifest <- fread(MANIFEST_FILE, skip = "IlmnID", select = c("IlmnID", "UCSC_RefGene_Name"))
manifest[, UCSC_RefGene_Name := fifelse(is.na(UCSC_RefGene_Name), "", UCSC_RefGene_Name)]

mirna_probes <- manifest[grepl("MIR", UCSC_RefGene_Name, ignore.case = TRUE), IlmnID]
cat("Found", length(mirna_probes), "miRNA-associated probes.\n\n")

# --- Step 2: Process each cancer ---
for (cancer in CANCERS) {
  cat("Processing:", cancer, "...")

  input_file  <- file.path(DATA_DIR, cancer, "M.csv")
  output_file <- file.path(OUTPUT_DIR, paste0("M_", cancer, "_miRNA_annotated.csv"))

  if (!file.exists(input_file)) { cat(" SKIP (M.csv not found)\n"); next }

  tryCatch({
    m_data    <- fread(input_file)
    probe_col <- colnames(m_data)[1]
    m_filtered <- m_data[get(probe_col) %in% mirna_probes]

    if (nrow(m_filtered) > 0) {
      fwrite(m_filtered, output_file)
      cat(" saved", nrow(m_filtered), "probes.\n")
    } else {
      cat(" no miRNA probes found.\n")
    }
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
  })
}

cat("\n=== miRNA Probe Extraction Complete ===\n")
