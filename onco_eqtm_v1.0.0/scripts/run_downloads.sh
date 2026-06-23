#!/bin/bash
# =============================================================================
# Run all Onco-eQTM data download scripts
#
# Downloads clinical, expression, methylation, miRNA, immune, pathway,
# and drug sensitivity data from public repositories.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOWNLOAD_DIR="${SCRIPT_DIR}/download"

cd "$PROJECT_ROOT"

echo "=== Onco-eQTM: Downloading all input data ==="
echo "Project root: ${PROJECT_ROOT}"
echo ""

bash "${DOWNLOAD_DIR}/01_download_clinical_data.sh"
bash "${DOWNLOAD_DIR}/02_download_gene_expression.sh"
bash "${DOWNLOAD_DIR}/03_download_methylation.sh"
bash "${DOWNLOAD_DIR}/04_download_mirna_expression.sh"
bash "${DOWNLOAD_DIR}/11_download_immune_signatures.sh"
bash "${DOWNLOAD_DIR}/12_download_paradigm_pathway.sh"
bash "${DOWNLOAD_DIR}/13_download_vaen_drug_data.sh"

echo ""
echo "All downloads complete."
