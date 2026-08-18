#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/06_annotations_frame_substrate.csv"
POINT="$ROOT/Data/Final/point_annotations.csv"
VME="$ROOT/Data/Final/vme_annotations.csv"
IMPL="$ROOT/Scripts/legacy_current/07_split_point_and_vme_CURRENT.sh"

mkdir -p "$ROOT/Data/Final"
for f in "$INPUT" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"

(cd "$TMP" && bash "$IMPL")
cp "$TMP/point_annotations.csv" "$POINT"
cp "$TMP/vme_annotations.csv" "$VME"

python3 - "$INPUT" "$POINT" "$VME" <<'PY'
import csv,sys
from pathlib import Path
paths=list(map(Path,sys.argv[1:]))
def load(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f); rows=list(r); return r.fieldnames or [],rows
h0,r0=load(paths[0]); hp,rp=load(paths[1]); hv,rv=load(paths[2])
if len(rp)+len(rv)!=len(r0): raise SystemExit("ERROR: split row counts do not sum to source")
if hp!=h0 or hv!=h0: raise SystemExit("ERROR: split headers differ from source")
if any((x.get("top_level") or "").strip().casefold()=="vme" for x in rp):
    raise SystemExit("ERROR: VME rows found in point_annotations.csv")
if any((x.get("top_level") or "").strip().casefold()!="vme" for x in rv):
    raise SystemExit("ERROR: non-VME rows found in vme_annotations.csv")
print(f"STEP 07 COMPLETE: point rows={len(rp):,}; VME rows={len(rv):,}")
print(f"  {paths[1]}")
print(f"  {paths[2]}")
PY
