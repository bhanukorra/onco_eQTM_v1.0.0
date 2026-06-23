#!/bin/bash
# =============================================================================
# Download and Process TCGA PARADIGM Pathway ssGSEA Data
#
# Loops through all 27 TCGA cancer types to download the corresponding 
# ssGSEA 1387 pathway scores (Z-normalized) from the UCSC Xena Hub.
# The downloaded text files are then compressed as expected by the analysis script.
# =============================================================================

set -euo pipefail

# Get directory where script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"

CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC \
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

BASE_URL="https://tcga.xenahubs.net/download/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level_Z"

echo "=== Downloading TCGA Pan-Cancer ssGSEA Pathway Data ==="

for CANCER in $CANCERS; do
    echo "Processing $CANCER..."
    cancer_dir="${DATA_DIR}/${CANCER}"
    mkdir -p "$cancer_dir"
    
    filename="${CANCER}_PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level_Z.txt"
    outfile="${cancer_dir}/${filename}"
    gzfile="${outfile}.gz"
    
    if [ -f "$gzfile" ]; then
        echo "  Parsed pathway file already exists at $gzfile. Skipping."
    else
        url="${BASE_URL}/${filename}"
        echo "  Downloading: $url"
        curl -L "$url" -o "$outfile"
        
        echo "  Compressing with gzip..."
        gzip -f "$outfile"
    fi
done

echo "Done."
