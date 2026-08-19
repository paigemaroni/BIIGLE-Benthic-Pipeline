#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/03_annotations_taxonomy.csv"
LOOKUP="$ROOT/Sheets/depth_temp_frameid.csv"
OUTPUT="$ROOT/Data/Intermediate/04_annotations_environment.csv"
IMPL="$ROOT/Scripts/internal/04_add_frame_metadata_impl.sh"

mkdir -p "$(dirname "$OUTPUT")"
for f in "$INPUT" "$LOOKUP" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"
cp "$LOOKUP" "$TMP/depth_temp_frameid.csv"

(cd "$TMP" && bash "$IMPL")
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"

python3 - "$INPUT" "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
a,b=map(Path,sys.argv[1:])
with a.open(encoding='utf-8-sig',newline='') as f: na=sum(1 for _ in csv.DictReader(f))
with b.open(encoding='utf-8-sig',newline='') as f:
    r=csv.DictReader(f); rows=list(r); h=r.fieldnames or []
if na!=len(rows): raise SystemExit("ERROR: row count changed")
for c in ("video_time","depth_m","temperature"):
    if c not in h: raise SystemExit(f"ERROR: output missing {c}")
matched=len({r["filename"] for r in rows if (r.get("depth_m") or "").strip() or (r.get("temperature") or "").strip()})
print(f"STEP 04 COMPLETE: {len(rows):,} rows; {matched:,} frames with depth/temp values -> {b}")
PY
