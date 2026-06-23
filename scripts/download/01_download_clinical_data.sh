#!/bin/bash
# Download TCGA Clinical Data from UCSC Xena Hub
# Downloads survival, age, sex and other clinical info for each cancer type.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${PROJECT_ROOT}/data"

CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

BASE_URL="https://tcga.xenahubs.net/download"

for CANCER in $CANCERS; do
    echo "Downloading clinical data for $CANCER..."
    mkdir -p "${OUTPUT_DIR}/${CANCER}"
    URL="${BASE_URL}/TCGA.${CANCER}.sampleMap/${CANCER}_clinicalMatrix"
    curl -L "$URL" -o "${OUTPUT_DIR}/${CANCER}/${CANCER}_clinicalMatrix"
done

echo "Done."
