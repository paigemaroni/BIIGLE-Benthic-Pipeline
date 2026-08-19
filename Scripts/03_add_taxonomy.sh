#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
INPUT="$ROOT/Data/Intermediate/02_annotations_top_level.csv"
LOOKUP="$ROOT/Sheets/taxon_list.csv"
OUTPUT="$ROOT/Data/Intermediate/03_annotations_taxonomy.csv"
QC="$ROOT/Analyses/00_Quality_Control"
IMPL="$ROOT/Scripts/internal/03_add_taxonomy_impl.sh"

mkdir -p "$(dirname "$OUTPUT")" "$QC"
for f in "$INPUT" "$LOOKUP" "$IMPL"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp "$INPUT" "$TMP/combined_biigle_annotations.csv"
cp "$LOOKUP" "$TMP/taxon_list.csv"

bash "$IMPL" "$TMP" "combined_biigle_annotations.csv" "taxon_list.csv"
cp "$TMP/combined_biigle_annotations.csv" "$OUTPUT"
[[ -f "$TMP/taxon_lookup_unmatched.csv" ]] && cp "$TMP/taxon_lookup_unmatched.csv" "$QC/03_taxon_lookup_unmatched.csv"
[[ -f "$TMP/taxon_lookup_partial_or_ambiguous.csv" ]] && cp "$TMP/taxon_lookup_partial_or_ambiguous.csv" "$QC/03_taxon_lookup_partial_or_ambiguous.csv"

# Split the raw unmatched report into:
#   - expected VME rows deferred to Step 05; and
#   - genuinely unmatched non-VME hierarchies requiring review.
if [[ -f "$QC/03_taxon_lookup_unmatched.csv" ]]; then
python3 - "$QC/03_taxon_lookup_unmatched.csv" "$QC" <<'PYQC'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
qc = Path(sys.argv[2])

with source.open(encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))

vme = []
non_vme = []

for row in rows:
    hierarchy = str(row.get("label_hierarchy", "") or "").strip()
    if hierarchy.casefold().startswith("vme"):
        vme.append(row)
    else:
        non_vme.append(row)

for path, values in [
    (qc / "03_taxon_lookup_deferred_vme.csv", vme),
    (qc / "03_taxon_lookup_unmatched_non_vme.csv", non_vme),
]:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["label_hierarchy", "row_count"]
        )
        writer.writeheader()
        writer.writerows(values)

print(
    "Step 03 unmatched hierarchy review: "
    f"{sum(int(x.get('row_count', 0) or 0) for x in vme):,} VME row(s) deferred to Step 05; "
    f"{sum(int(x.get('row_count', 0) or 0) for x in non_vme):,} non-VME row(s) remain unmatched."
)
PYQC
fi

python3 - "$INPUT" "$OUTPUT" <<'PY'
import csv, sys
from pathlib import Path
a,b=map(Path,sys.argv[1:])
with a.open(encoding='utf-8-sig',newline='') as f:
    ra=list(csv.DictReader(f))
with b.open(encoding='utf-8-sig',newline='') as f:
    r=csv.DictReader(f); hb=r.fieldnames or []; rb=list(r)
required=["cpc_codes","kingdom","phylum","class","order","family","taxonomic_resolution","morphology","common_id_short","common_id_mid","common_id_full","releif","substrate","type","size"]
if len(ra)!=len(rb): raise SystemExit("ERROR: row count changed")
missing=[x for x in required if x not in hb]
if missing: raise SystemExit("ERROR: taxonomy output missing "+", ".join(missing))
print(f"STEP 03 COMPLETE: {len(rb):,} rows -> {b}")
PY
