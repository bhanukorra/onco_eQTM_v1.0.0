#!/bin/bash
# =============================================================================
# Run all Onco-eQTM analysis scripts
#
# Executes core eQTM discovery and all downstream association analyses.
# Run scripts/run_downloads.sh first to fetch input data.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ANALYSIS_DIR="${SCRIPT_DIR}/analysis"

cd "$PROJECT_ROOT"

echo "=== Onco-eQTM: Running all analyses ==="
echo "Project root: ${PROJECT_ROOT}"
echo ""

echo "--- Core eQTM discovery ---"
bash "${ANALYSIS_DIR}/05_run_eqtm_pipeline.sh"

echo ""
echo "--- Downstream analyses ---"
Rscript "${ANALYSIS_DIR}/06_meth_drug_correlation.R"
Rscript "${ANALYSIS_DIR}/07_immune_score_analysis.R"
Rscript "${ANALYSIS_DIR}/08_mirna_probe_extraction.R"
Rscript "${ANALYSIS_DIR}/09_mirna_eqtl_analysis.R"
Rscript "${ANALYSIS_DIR}/10_paradigm_analysis.R"

echo ""
echo "All analyses complete."
