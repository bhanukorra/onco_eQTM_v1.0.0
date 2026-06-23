#!/bin/bash
# =============================================================================
# Download and Process VAEN Imputed Drug Response Data
#
# This script downloads the imputed drug response values (ln(IC50)) for TCGA 
# cancer samples from the VAEN bioinformatics repository and splits the data 
# into cancer-specific files (e.g. {CANCER}_drug_data.csv).
#
# Options:
#   -d <dataset>  Specify dataset to download: "GDSC" (default, 251 drugs) 
#                 or "CCLE" (24 drugs).
# =============================================================================

set -euo pipefail

# Default option
DATASET="GDSC"

# Parse arguments
while getopts "d:" opt; do
  case $opt in
    d)
      DATASET=$(echo "$OPTARG" | tr '[:lower:]' '[:upper:]')
      ;;
    *)
      echo "Usage: $0 [-d GDSC|CCLE]" >&2
      exit 1
      ;;
  esac
done

if [ "$DATASET" != "GDSC" ] && [ "$DATASET" != "CCLE" ]; then
    echo "ERROR: Invalid dataset specified. Must be 'GDSC' or 'CCLE'." >&2
    exit 1
fi

# Get directory where script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"

echo "=== Downloading VAEN Imputed Drug Response Data ($DATASET) ==="
mkdir -p "$DATA_DIR"

if [ "$DATASET" = "GDSC" ]; then
    URL="https://bioinfo.uth.edu/VAEN/result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt"
    OUTPUT_FILE="${DATA_DIR}/VAEN_GDSC.A.pred_TCGA.txt"
else
    URL="https://bioinfo.uth.edu/VAEN/result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt"
    OUTPUT_FILE="${DATA_DIR}/VAEN_CCLE.A.pred_TCGA.txt"
fi

if [ -f "$OUTPUT_FILE" ]; then
    echo "VAEN drug response file already exists at $OUTPUT_FILE. Skipping download."
else
    echo "Downloading from: $URL"
    curl -L "$URL" -o "$OUTPUT_FILE"
fi

echo "Splitting drug response predictions by cancer type..."
python3 "${SCRIPT_DIR}/split_vaen_drug_data.py"

echo "Done."
