#!/bin/bash
# Download TCGA Gene Expression Data (HiSeqV2) from UCSC Xena Hub
# RNA-Seq gene expression normalized as log2(count+1).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${PROJECT_ROOT}/data"

CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

BASE_URL="https://tcga.xenahubs.net/download"

for CANCER in $CANCERS; do
    echo "Downloading gene expression for $CANCER..."
    mkdir -p "${OUTPUT_DIR}/${CANCER}"
    URL="${BASE_URL}/TCGA.${CANCER}.sampleMap%2FHiSeqV2.gz"
    wget -c "$URL" -O "${OUTPUT_DIR}/${CANCER}/${CANCER}_HiSeqV2.gz"
done

echo "Done."
