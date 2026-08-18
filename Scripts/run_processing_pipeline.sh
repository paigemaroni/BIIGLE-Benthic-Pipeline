#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"

echo "BIIGLE processing pipeline"
echo "=========================="
echo "Project: $ROOT"
echo

for script in \
    00_validate_inputs.sh \
    01_add_lat_long.sh \
    02_add_top_level.sh \
    03_add_taxonomy.sh \
    04_add_frame_metadata.sh \
    05_add_vme_taxonomy.sh \
    06_add_frame_substrate_class.sh \
    07_split_point_and_vme_data.sh
do
    echo
    echo ">>> $script"
    bash "$ROOT/Scripts/$script" "$ROOT"
done

python3 - "$ROOT" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
files = [
    ("01", "latitude/longitude", root / "Data/Intermediate/01_annotations_latlong.csv"),
    ("02", "top-level class", root / "Data/Intermediate/02_annotations_top_level.csv"),
    ("03", "taxonomy/abiotic enrichment", root / "Data/Intermediate/03_annotations_taxonomy.csv"),
    ("04", "frame environmental metadata", root / "Data/Intermediate/04_annotations_environment.csv"),
    ("05", "VME taxonomy", root / "Data/Intermediate/05_annotations_vme_taxonomy.csv"),
    ("06", "frame substrate", root / "Data/Intermediate/06_annotations_frame_substrate.csv"),
    ("07a", "final point data", root / "Data/Final/point_annotations.csv"),
    ("07b", "final VME data", root / "Data/Final/vme_annotations.csv"),
]
out = root / "Analyses/00_Quality_Control/07_processing_pipeline_manifest.csv"
out.parent.mkdir(parents=True, exist_ok=True)

rows = []
for stage, description, path in files:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        n_rows = sum(1 for _ in reader)
    rows.append([stage, description, str(path.relative_to(root)), n_rows, len(header)])

with out.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["stage", "description", "output_file", "rows", "columns"])
    writer.writerows(rows)

print(f"Processing manifest: {out}")
PY

echo
echo "PROCESSING PIPELINE COMPLETE"
echo "Final data:"
echo "  $ROOT/Data/Final/point_annotations.csv"
echo "  $ROOT/Data/Final/vme_annotations.csv"
