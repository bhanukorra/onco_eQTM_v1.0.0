#!/usr/bin/env Rscript
# =============================================================================
# Pan-Cancer miRNA cis-eQTL Analysis
#
# Runs MatrixEQTL to find cis-eQTLs between methylation probes (CpG sites)
# and miRNA expression levels across 27 cancer types.
# It also annotates mature miRNA targets using a GFF3 annotation file.
# =============================================================================

if (!require("MatrixEQTL")) install.packages("MatrixEQTL")
if (!require("data.table")) install.packages("data.table")
library("MatrixEQTL")
library("data.table")

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
DATA_DIR                <- "./data"

# miRNA Annotations GFF3 file (e.g. hsa.gff3 from miRBase)
GFF3_ANNOTATION_PATH    <- "./data/hsa.gff3"

# Location files containing genomic coordinates
GENE_LOCATION_FILE      <- "./data/geneloc_miRNA.txt"
SNPS_LOCATION_FILE      <- "./data/snpsloc.txt"

# Where to save the output files
OUTPUT_BASE_DIR         <- "./results/mirna_eqtl"

# =============================================================================

# Parameters
fdr_thresholds <- c(0.05) 
useModel = modelLINEAR
pvOutputThreshold.cis = 1
pvOutputThreshold.trans = 0
cisDist = 1e6
errorCovariance = numeric()

# 27 TCGA cancer types
CANCERS <- c(
  "ACC", "BLCA", "BRCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UVM"
)

# =============================================================================
# STEP 1: Load Annotations
# =============================================================================
message("--------------------------------------------------")
message("Step 1: Loading miRNA Annotations...")

if (file.exists(GFF3_ANNOTATION_PATH)) {
  gff_data <- fread(GFF3_ANNOTATION_PATH, header = FALSE, sep = "\t", skip = "#", select = c(3, 9))
  gff_miRNA <- gff_data[V3 == "miRNA"]
  gff_miRNA[, gene := sub(".*ID=(MIMAT[0-9]+).*", "\\1", V9)]
  gff_miRNA[, miRNA_Name := sub(".*Name=([^;]+).*", "\\1", V9)]
  
  miRNA_map <- unique(gff_miRNA[, .(gene, miRNA_Name)], by = "gene")
  message(paste("  Loaded annotations for", nrow(miRNA_map), "unique miRNAs."))
} else {
  warning(paste("  WARNING: GFF3 file not found at", GFF3_ANNOTATION_PATH))
  message("  Results will NOT be annotated with names.")
  miRNA_map <- NULL
}

# =============================================================================
# STEP 2: Pan-Cancer Analysis Loop
# =============================================================================
for (cancer in CANCERS) {
  
  message(paste("\n========================================================"))
  message(paste("Processing Cancer:", cancer))
  
  # Define Input / Output Paths
  cancer_dir     <- file.path(DATA_DIR, cancer)
  cis_output_dir <- file.path(OUTPUT_BASE_DIR, cancer)
  
  SNP_file_name        <- file.path(cancer_dir, "M.csv")
  expression_file_name <- file.path(cancer_dir, "G.csv") 
  covariates_file_name <- file.path(cancer_dir, "C.csv")
  
  # Validate Input Files
  if (!all(file.exists(SNP_file_name, expression_file_name, covariates_file_name))) {
    message(paste("  SKIP: Missing M.csv, G.csv, or C.csv for", cancer))
    next
  }
  
  if (!dir.exists(cis_output_dir)) dir.create(cis_output_dir, recursive = TRUE)
  
  tryCatch({
    # --- Align Samples ---
    message("  Synchronizing samples...")
    snp_data  <- fread(SNP_file_name, header = TRUE)
    expr_data <- fread(expression_file_name, header = TRUE)
    cov_data  <- fread(covariates_file_name, header = TRUE)
    
    snp_samples  <- colnames(snp_data)[-1]
    expr_samples <- colnames(expr_data)[-1]
    cov_samples  <- cov_data[[1]]
    
    common_samples <- intersect(intersect(snp_samples, expr_samples), cov_samples)
    if (length(common_samples) == 0) stop("No matching samples found.")
    message(paste("    Matched", length(common_samples), "samples."))
    
    # Align matrices
    cols_snp <- c(colnames(snp_data)[1], common_samples)
    snp_data_aligned <- snp_data[, ..cols_snp]
    
    cols_expr <- c(colnames(expr_data)[1], common_samples)
    expr_data_aligned <- expr_data[, ..cols_expr]
    
    cov_data_aligned <- cov_data[cov_data[[1]] %in% common_samples, ]
    cov_data_aligned <- cov_data_aligned[match(common_samples, cov_data_aligned[[1]]), ]
    
    cov_mat <- t(as.matrix(cov_data_aligned[, -1, with = FALSE]))
    colnames(cov_mat) <- cov_data_aligned[[1]]
    
    # Write Temp Files
    temp_snp_file  <- file.path(cis_output_dir, "M_aligned_temp.csv")
    temp_expr_file <- file.path(cis_output_dir, "G_aligned_temp.csv")
    temp_cov_file  <- file.path(cis_output_dir, "C_aligned_temp.csv")
    
    fwrite(snp_data_aligned, temp_snp_file)
    fwrite(expr_data_aligned, temp_expr_file)
    write.csv(cov_mat, temp_cov_file, quote = FALSE)
    
    # --- Run Matrix eQTL ---
    message("  Running MatrixEQTL...")
    snps = SlicedData$new(); snps$fileDelimiter = ","; snps$fileOmitCharacters = "NA"; snps$fileSkipRows = 1; snps$fileSkipColumns = 1; snps$fileSliceSize = 2000
    snps$LoadFile(temp_snp_file)
    
    gene = SlicedData$new(); gene$fileDelimiter = ","; gene$fileOmitCharacters = "NA"; gene$fileSkipRows = 1; gene$fileSkipColumns = 1; gene$fileSliceSize = 2000
    gene$LoadFile(temp_expr_file)
    
    cvrt = SlicedData$new(); cvrt$fileDelimiter = ","; cvrt$fileOmitCharacters = "NA"; cvrt$fileSkipRows = 1; cvrt$fileSkipColumns = 1
    cvrt$LoadFile(temp_cov_file)
    
    snpspos = fread(SNPS_LOCATION_FILE)
    setnames(snpspos, c("snpid", "chr", "pos"))
    
    genepos = fread(GENE_LOCATION_FILE)
    setnames(genepos, c("geneid", "chr", "s1", "s2"))
    
    output_file_name_cis_raw = file.path(cis_output_dir, "eQTL_results_cis_RAW.csv")
    
    me = Matrix_eQTL_main(
      snps = snps, gene = gene, cvrt = cvrt,
      output_file_name.cis = output_file_name_cis_raw,
      output_file_name = NULL,
      pvOutputThreshold.cis = pvOutputThreshold.cis,
      pvOutputThreshold = pvOutputThreshold.trans,
      useModel = useModel, errorCovariance = errorCovariance,
      snpspos = snpspos, genepos = genepos, cisDist = cisDist,
      pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE, verbose = FALSE
    )
    
    # --- Annotate and Save ---
    message("  Processing results...")
    if (file.exists(output_file_name_cis_raw) && file.size(output_file_name_cis_raw) > 0) {
      res_cis <- fread(output_file_name_cis_raw)
      
      if (nrow(res_cis) > 0) {
        if (!is.null(miRNA_map)) {
          message("    Annotating with miRNA names...")
          res_cis <- merge(res_cis, miRNA_map, by = "gene", all.x = TRUE)
          if ("miRNA_Name" %in% names(res_cis)) {
            col_order     <- c("SNP", "gene", "miRNA_Name", "beta", "t-stat", "p-value", "FDR")
            existing_cols <- intersect(col_order, names(res_cis))
            res_cis       <- res_cis[, ..existing_cols]
          }
        }
        
        for (fdr in fdr_thresholds) {
          res_filtered <- res_cis[FDR < fdr]
          fdr_file     <- file.path(cis_output_dir, paste0(cancer, "_cis_mirna_FDR_", fdr, "_annotated.csv"))
          
          if (nrow(res_filtered) > 0) {
            res_filtered <- res_filtered[order(FDR)]
            fwrite(res_filtered, fdr_file)
            message(paste("    -> SAVED:", nrow(res_filtered), "annotated eQTLs at FDR <", fdr))
            
            # Save Max P-value info
            max_p      <- max(res_filtered$`p-value`)
            max_p_file <- file.path(cis_output_dir, paste0("Max_PVal_FDR_", fdr, ".txt"))
            writeLines(paste("Max P-value:", max_p), max_p_file)
          } else {
            message(paste("    -> No eQTLs found at FDR <", fdr))
          }
        }
      } else {
        message("    MatrixEQTL ran but found 0 results.")
      }
    }
    
    # Cleanup Temp Files
    unlink(c(temp_snp_file, temp_expr_file, temp_cov_file))
    
  }, error = function(e) {
    message(paste("  ERROR processing", cancer, ":", e$message))
  })
}

message("\n--------------------------------------------------")
message("miRNA cis-eQTL Analysis Complete.")
