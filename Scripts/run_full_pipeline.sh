#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"

bash "$ROOT/Scripts/run_processing_pipeline.sh" "$ROOT"
bash "$ROOT/Scripts/run_analysis_pipeline.sh" "$ROOT"
