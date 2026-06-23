#!/bin/bash
# Download TCGA miRNA Expression Data (miRNA HiSeq) from UCSC Xena Hub
# miRNA gene-level expression values.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${PROJECT_ROOT}/data"

CANCERS="ACC BLCA BRCA CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC
         LUAD LUSC MESO PAAD PCPG READ SARC SKCM STAD TGCT THCA THYM UVM"

BASE_URL="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"

for CANCER in $CANCERS; do
    echo "Downloading miRNA expression for $CANCER..."
    mkdir -p "${OUTPUT_DIR}/${CANCER}"
    URL="${BASE_URL}/TCGA.${CANCER}.sampleMap%2FmiRNA_HiSeq_gene.gz"
    wget -c "$URL" -O "${OUTPUT_DIR}/${CANCER}/${CANCER}_miRNA_HiSeq_gene.gz"
done

echo "Done."
