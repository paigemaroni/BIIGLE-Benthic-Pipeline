#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/05_annotations_vme_taxonomy.csv"
OUTPUT="$ROOT/Data/Intermediate/06_annotations_frame_substrate.csv"
ANALYSIS_DIR="$ROOT/Analyses/01_Frame_Composition"
IMPL="$ROOT/Scripts/legacy_current/06_add_frame_substrate_class_CURRENT.sh"

mkdir -p "$(dirname "$OUTPUT")" "$ANALYSIS_DIR"
for f in "$INPUT" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"

(cd "$TMP" && bash "$IMPL")
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"
[[ -f "$TMP/frame_substrate_summary.csv" ]] && cp "$TMP/frame_substrate_summary.csv" "$ANALYSIS_DIR/06_frame_substrate_summary.csv"

python3 - "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
p=Path(sys.argv[1])
with p.open(encoding='utf-8-sig',newline='') as f:
    r=csv.DictReader(f); rows=list(r); h=r.fieldnames or []
for c in ("frame_substrate_class","frame_relief_class","frame_substrate_n","frame_substrate_dominance_pct"):
    if c not in h: raise SystemExit(f"ERROR: output missing {c}")
frames={}
for row in rows:
    frames.setdefault(row["filename"], (row.get("frame_substrate_class",""), row.get("frame_relief_class","")))
classified=sum(1 for s,r in frames.values() if str(s).strip())
print(f"STEP 06 COMPLETE: {classified:,}/{len(frames):,} frames with substrate class -> {p}")
PY
