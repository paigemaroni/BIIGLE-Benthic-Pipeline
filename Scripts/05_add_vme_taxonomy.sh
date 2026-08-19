#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/04_annotations_environment.csv"
LOOKUP="$ROOT/Sheets/vme_taxon_list.csv"
OUTPUT="$ROOT/Data/Intermediate/05_annotations_vme_taxonomy.csv"
QC="$ROOT/Analyses/00_Quality_Control"
IMPL="$ROOT/Scripts/internal/05_add_vme_taxonomy_impl.sh"

mkdir -p "$(dirname "$OUTPUT")" "$QC"
for f in "$INPUT" "$LOOKUP" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"
cp "$LOOKUP" "$TMP/vme_taxon_list.csv"

(cd "$TMP" && bash "$IMPL")
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"
[[ -f "$TMP/vme_taxonomy_unmatched.csv" ]] && cp "$TMP/vme_taxonomy_unmatched.csv" "$QC/05_vme_taxonomy_unmatched.csv"
[[ -f "$TMP/vme_taxonomy_ambiguous.csv" ]] && cp "$TMP/vme_taxonomy_ambiguous.csv" "$QC/05_vme_taxonomy_ambiguous.csv"

python3 - "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
p=Path(sys.argv[1])
with p.open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f))
v=[r for r in rows if (r.get("top_level") or "").strip().casefold()=="vme"]
print(f"STEP 05 COMPLETE: {len(v):,} VME rows processed; see 05_vme_taxonomy_unmatched.csv and 05_vme_taxonomy_ambiguous.csv for unresolved records -> {p}")
PY
