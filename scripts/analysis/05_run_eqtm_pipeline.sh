#!/bin/bash
# =============================================================================
# Pan-Cancer eQTM Analysis Pipeline
#
#
# Output for each cancer:
#   Results/final_CANCER.csv    → spreadsheet of significant CpG-gene pairs
#   Results/final_CANCER.sqlite → same data as a searchable database
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# =============================================================================
# SETTINGS — change these to match your system
# =============================================================================

# Where your cancer data folders are (each should have M.csv, G.csv, C.csv)
DATA_DIR="${PROJECT_ROOT}/data"

# Where Torch-eCpG is installed
TECPG_DIR="${PROJECT_ROOT}/torch-ecpg"

# =============================================================================
# Everything below runs automatically
# =============================================================================

TECPG_ANALYSIS="${TECPG_DIR}/Analysis"
GENE_BED="${TECPG_ANALYSIS}/annot/G.bed6"
METH_BED="${TECPG_ANALYSIS}/annot/M.bed6"
REGION_PY="${TECPG_DIR}/demo/assignRegionToEcpg.py"
LOG="${DATA_DIR}/pipeline_progress.log"

export PYTHONPATH="${TECPG_DIR}:${PYTHONPATH:-}"

# 27 TCGA cancer types
#   ACC(79) BLCA(407) BRCA(783) CHOL(36) COAD(278) DLBC(48) ESCA(184) GBM(51)
#   HNSC(520) KICH(66) KIRC(318) KIRP(274) LGG(516) LIHC(371) LUAD(453)
#   LUSC(371) MESO(87) PAAD(178) PCPG(179) READ(92) SARC(259) SKCM(103)
#   STAD(372) TGCT(150) THCA(505) THYM(120) UVM(80)
CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

[ ! -d "$DATA_DIR" ]   && echo "ERROR: Data folder not found: $DATA_DIR" && exit 1
[ ! -d "$TECPG_DIR" ]  && echo "ERROR: Torch-eCpG not found: $TECPG_DIR" && exit 1

echo "Pipeline started at $(date)" | tee "$LOG"

for CANCER in $CANCERS; do
    DIR="${DATA_DIR}/${CANCER}"
    [[ ! -d "$DIR" ]] && continue
    [[ ! -f "$DIR/M.csv" || ! -f "$DIR/G.csv" || ! -f "$DIR/C.csv" ]] && \
        echo "[SKIP] $CANCER: missing input files" | tee -a "$LOG" && continue

    echo ""
    echo ">>> Processing: $CANCER"

    RESULTS="${DIR}/Results"
    TORCH_OUT="${RESULTS}/torch_output"
    MATRIX_OUT="${RESULTS}/matrix_eqtl_output"
    mkdir -p "$RESULTS" "$TORCH_OUT" "$MATRIX_OUT"

    # ---- Step 0: Clean covariates (remove constant columns) ----
    echo "  [0] Cleaning covariates..."
    Rscript -e "
        library(data.table)
        df <- fread('${DIR}/C.csv', header = TRUE)
        cols <- names(df)[-1]
        bad <- cols[sapply(df[, ..cols], function(x) length(unique(x)) <= 1)]
        if (length(bad) > 0) df <- df[, !bad, with = FALSE]
        fwrite(df, '${DIR}/C_filtered.csv')
    "

    # ---- Step 1: Run Torch-eCpG ----
    echo "  [1] Running Torch-eCpG..."
    TMP="${DIR}/.tecpg_input" && mkdir -p "$TMP"
    ln -sf "$DIR/M.csv" "$TMP/M.csv"
    ln -sf "$DIR/G.csv" "$TMP/G.csv"
    ln -sf "$DIR/C_filtered.csv" "$TMP/C.csv"

    python3 -m tecpg \
        --root-path "$TECPG_ANALYSIS" --input-dir "$TMP" --output-dir "$TORCH_OUT" \
        --meth-file M.csv --gene-file G.csv --covar-file C.csv \
        --meth-annot "$METH_BED" --gene-annot "$GENE_BED" \
        --cpu-threads 16 run mlr --cis -p 1 -g 10000 -m 10000
    rm -rf "$TMP"

    # Merge chunks
    cd "$TORCH_OUT"
    shopt -s nullglob; files=(*-*.csv)
    [ ${#files[@]} -gt 0 ] && \
        (head -n 1 "${files[0]}"; tail -q -n +2 "${files[@]}" | sort -V) > merged_tecpg_output.csv && \
        rm -f *-*.csv

    # FDR correction (keep FDR <= 0.05)
    Rscript -e "
        library(data.table)
        df <- fread('merged_tecpg_output.csv')
        df\$fdr <- p.adjust(df\$mt_p, method = 'BH')
        cutoff <- max(df\$mt_p[df\$fdr <= 0.05], na.rm = TRUE)
        if (is.infinite(cutoff)) cutoff <- 0
        fwrite(df[df\$mt_p <= cutoff, ], 'significant_fdr_0.05.csv')
    "

    # Annotate genomic regions
    python3 "$REGION_PY" -d "significant_fdr_0.05.csv" -g "$GENE_BED" -m "$METH_BED" \
        -o "significant_fdr_0.05_annotated.csv"

    # ---- Step 2: Run Matrix eQTL ----
    echo "  [2] Running Matrix eQTL..."
    Rscript -e "
        library(MatrixEQTL); library(data.table)
        snps = SlicedData\$new(); snps\$fileDelimiter = ','; snps\$LoadFile('${DIR}/M.csv')
        gene = SlicedData\$new(); gene\$fileDelimiter = ','; gene\$LoadFile('${DIR}/G.csv')
        cvrt_raw <- fread('${DIR}/C_filtered.csv')
        ids <- cvrt_raw[[1]]
        trans_dt <- as.data.table(transpose(cvrt_raw[, -1, with=FALSE]))
        setnames(trans_dt, names(trans_dt), ids)
        trans_dt[, id := names(cvrt_raw)[-1]]; setcolorder(trans_dt, 'id')
        fwrite(trans_dt, '${MATRIX_OUT}/transformed_covariates.csv', row.names=FALSE)
        cvrt = SlicedData\$new(); cvrt\$fileDelimiter = ','; cvrt\$LoadFile('${MATRIX_OUT}/transformed_covariates.csv')
        write(gene\$nCols(), '${MATRIX_OUT}/n_samples.txt')
        write(cvrt\$nRows(), '${MATRIX_OUT}/n_covariates.txt')
        snpspos = fread('${METH_BED}')[, c(4,1,2)]; setnames(snpspos, c('snpid','chr','pos'))
        genepos = fread('${GENE_BED}')[, c(4,1,2,3)]; setnames(genepos, c('geneid','chr','s1','s2'))
        me = Matrix_eQTL_main(snps=snps, gene=gene, cvrt=cvrt, output_file_name=NULL, pvOutputThreshold=0, output_file_name.cis='${MATRIX_OUT}/eQTL_results_cis_RAW.csv', pvOutputThreshold.cis=1, useModel=modelLINEAR, snpspos=snpspos, genepos=genepos, cisDist=1e6, verbose=TRUE, pvalue.hist=FALSE)
        res <- me\$cis\$eqtls
        res\$FDR <- p.adjust(res\$pvalue, method='BH')
        cutoff <- max(res\$pvalue[res\$FDR <= 0.05], na.rm=TRUE)
        if (is.infinite(cutoff)) cutoff <- 0
        fwrite(res[res\$pvalue <= cutoff, ], '${MATRIX_OUT}/eQTL_results_cis_FINAL_FDR0.05.csv')
    "

    # ---- Step 3: Merge both results → final consensus ----
    echo "  [3] Building consensus output..."
    FINAL_CSV="${RESULTS}/final_${CANCER}.csv"
    Rscript -e "
        library(data.table); library(dplyr)
        torch <- fread('significant_fdr_0.05.csv')
        annot <- fread('significant_fdr_0.05_annotated.csv')
        torch_combined <- merge(annot, torch, by = c('mt_id','gt_id'), all.x = TRUE)
        matrix <- fread('${MATRIX_OUT}/eQTL_results_cis_FINAL_FDR0.05.csv')
        n_s <- as.numeric(readLines('${MATRIX_OUT}/n_samples.txt'))
        n_c <- as.numeric(readLines('${MATRIX_OUT}/n_covariates.txt'))
        df_val <- n_s - 2 - n_c
        torch_ren <- torch_combined %>% rename(SNP = mt_id, gene = gt_id)
        common <- inner_join(matrix, torch_ren, by = c('snps'='SNP','gene'='gene'))
        common\$correlation <- common\$mt_t / sqrt(df_val + common\$mt_t^2)
        strong <- common[abs(common\$correlation) > 0.3, ]
        final <- strong %>%
          select(mt_id=snps, gt_id=gene, mt_chrom, mt_chromStart, mt_strand, gt_chrom, gt_chromStart, gt_strand, region, mt_est, mt_err, mt_t, mt_p, beta, 't-stat'=statistic, 'p-value'=pvalue, FDR=FDR_manual, correlation) %>%
          mutate(rp = case_when(region=='CIS'~1, region=='PROMOTER'~2, region=='DISTAL'~3, TRUE~4)) %>%
          arrange(mt_id, gt_id, rp) %>% group_by(mt_id, gt_id) %>% slice(1) %>% ungroup() %>% select(-rp)
        fwrite(final, '${FINAL_CSV}')
    "
    echo "[DONE] $CANCER" | tee -a "$LOG"
done

echo ""
echo "All cancers processed. Results in each cancer's Results/ folder."
