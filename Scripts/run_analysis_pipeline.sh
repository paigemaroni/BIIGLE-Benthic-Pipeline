#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
cd "$ROOT"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript was not found." >&2
  echo "Install R, then run: Rscript Scripts/install_R_packages.R" >&2
  exit 1
fi

echo "BIIGLE ecological analysis pipeline"
echo "==================================="
echo "Project: $ROOT"
echo

for script in \
  20_analysis_qc.R \
  21_analysis_alpha_diversity.R \
  22_analysis_nmds.R \
  23_analysis_permanova.R \
  24_analysis_cap_dbrda.R \
  25_analysis_environment.R \
  26_analysis_univariate.R \
  27_analysis_taxon_responses.R \
  28_analysis_vme.R
do
  echo
  echo ">>> $script"
  Rscript "Scripts/$script"
done

echo
echo "ANALYSIS PIPELINE COMPLETE"
echo "Results: $ROOT/Analyses"
echo "Figures: $ROOT/Figures"
