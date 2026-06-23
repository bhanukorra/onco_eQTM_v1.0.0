#!/bin/bash
# Download TCGA DNA Methylation Data (HumanMethylation450) from UCSC Xena Hub
# Illumina 450K methylation beta values.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${PROJECT_ROOT}/data"

CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

BASE_URL="https://tcga.xenahubs.net/download"

for CANCER in $CANCERS; do
    echo "Downloading methylation data for $CANCER..."
    mkdir -p "${OUTPUT_DIR}/${CANCER}"
    URL="${BASE_URL}/TCGA.${CANCER}.sampleMap%2FHumanMethylation450.gz"
    wget -c "$URL" -O "${OUTPUT_DIR}/${CANCER}/${CANCER}_HumanMethylation450.gz"
done

echo "Done."
