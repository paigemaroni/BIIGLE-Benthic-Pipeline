#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/01_annotations_latlong.csv"
OUTPUT="$ROOT/Data/Intermediate/02_annotations_top_level.csv"
IMPL="$ROOT/Scripts/internal/02_add_top_level_impl.sh"

mkdir -p "$(dirname "$OUTPUT")"
for f in "$INPUT" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"

bash "$IMPL" "$TMP" "combined_biigle_annotations.csv"
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"

python3 - "$INPUT" "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
a,b=map(Path,sys.argv[1:])
def load(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f); rows=list(r); return r.fieldnames or [], rows
ha,ra=load(a); hb,rb=load(b)
if len(ra)!=len(rb): raise SystemExit("ERROR: row count changed")
if "top_level" not in hb: raise SystemExit("ERROR: top_level missing")
blank=sum(1 for x in rb if not (x.get("top_level") or "").strip())
print(f"STEP 02 COMPLETE: {len(rb):,} rows; blank top_level={blank:,} -> {b}")
PY
