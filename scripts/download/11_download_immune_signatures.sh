#!/bin/bash
# =============================================================================
# Download and Process TCGA Pan-Cancer Immune Infiltration Signatures
#
# This script downloads the 68 immune signatures for 10,852 whitelist samples
# from the Pan-Cancer Atlas Hub and splits the data into cancer-specific files.
# =============================================================================

set -euo pipefail

# Get directory where script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"

echo "=== Downloading TCGA Pan-Cancer 68 Immune Signatures ==="
mkdir -p "$DATA_DIR"

URL="https://pancanatlas.xenahubs.net/download/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz"
OUTPUT_FILE="${DATA_DIR}/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz"

if [ -f "$OUTPUT_FILE" ]; then
    echo "Immune signatures file already exists at $OUTPUT_FILE. Skipping download."
else
    echo "Downloading from: $URL"
    curl -L "$URL" -o "$OUTPUT_FILE"
fi

echo "Splitting immune signatures by cancer type..."
python3 "${SCRIPT_DIR}/split_immune_signatures.py"

echo "Done."
