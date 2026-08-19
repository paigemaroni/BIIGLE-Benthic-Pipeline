#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Sheets/combined_biigle_annotations.csv"
LOOKUP="$ROOT/Sheets/gl_latlong.csv"
OUTPUT="$ROOT/Data/Intermediate/01_annotations_latlong.csv"
IMPL="$ROOT/Scripts/internal/01_add_lat_long_impl.sh"

mkdir -p "$(dirname "$OUTPUT")"
for f in "$INPUT" "$LOOKUP" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"
cp "$LOOKUP" "$TMP/gl_latlong.csv"

bash "$IMPL" "$TMP" "combined_biigle_annotations.csv" "gl_latlong.csv"
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"

python3 - "$INPUT" "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
a,b=map(Path,sys.argv[1:])
def info(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        r=csv.reader(f); h=next(r); n=sum(1 for _ in r)
    return h,n
ha,na=info(a); hb,nb=info(b)
if na != nb: raise SystemExit(f"ERROR: row count changed: {na} -> {nb}")
for c in ("lat","long"):
    if c not in hb: raise SystemExit(f"ERROR: output missing {c}")
print(f"STEP 01 COMPLETE: {nb:,} rows -> {b}")
PY
