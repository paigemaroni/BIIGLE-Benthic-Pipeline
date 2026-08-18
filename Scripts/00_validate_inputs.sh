#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"

python3 - "$ROOT" <<'PY2'
import csv
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1]).resolve()
sheets = root / "Sheets"
out = root / "Analyses" / "00_Quality_Control"
out.mkdir(parents=True, exist_ok=True)

TARGET_POINTS = 100

required = {
    "combined_biigle_annotations.csv": [
        "dive_id", "label_name", "label_hierarchy", "filename",
        "shape_name", "annotation_id"
    ],
    "gl_latlong.csv": ["FullID", "Lat", "Long"],
    "depth_temp_frameid.csv": ["FrameID", "video_time", "depth_m", "temperature"],
    "taxon_list.csv": [
        "Layer_01", "cpc_codes", "kingdom", "phylum", "class",
        "taxonomic_resolution", "common_id_short", "common_id_mid",
        "common_id_full", "releif", "substrate", "type", "size"
    ],
    "vme_taxon_list.csv": [
        "Layer_01", "cpc_codes", "kingdom", "phylum", "class",
        "taxonomic_resolution", "common_id_short", "common_id_mid",
        "common_id_full"
    ],
}

def clean(value):
    return str(value or "").replace("\ufeff", "").strip()

def n(value):
    return clean(value).casefold()

def norm_id(value):
    value = clean(value).replace("\\", "/").split("/")[-1]
    value = re.sub(r"\.csv$", "", value, flags=re.I)
    m = re.search(r"20\d\d.*$", value)
    if m:
        value = m.group(0)
    value = value.replace("–", "-").replace("—", "-").replace("_", "-")
    value = re.sub(r"\s+", "", value)
    return value.casefold()

def top_level(hierarchy):
    hierarchy = clean(hierarchy)
    return hierarchy.split(">", 1)[0].strip() if hierarchy else ""

report = []
fatal = False

# -------------------------------------------------------------------------
# Required files and columns
# -------------------------------------------------------------------------
for filename, columns in required.items():
    path = sheets / filename
    if not path.exists():
        report.append((filename, "ERROR", "missing file"))
        fatal = True
        continue

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        row_count = sum(1 for _ in reader)

    normalised_header = {n(x) for x in header}
    missing = [col for col in columns if n(col) not in normalised_header]

    if missing:
        fatal = True

    report.append((
        filename,
        "PASS" if not missing else "ERROR",
        f'{row_count:,} rows; missing required columns: '
        f'{", ".join(missing) if missing else "none"}'
    ))

if fatal:
    with (out / "00_input_validation_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["check", "status", "details"])
        writer.writerows(report)

    print("BIIGLE input validation")
    print("=======================")
    for check, status, details in report:
        print(f"{status:7}  {check}: {details}")
    raise SystemExit(1)

# -------------------------------------------------------------------------
# Load principal files
# -------------------------------------------------------------------------
with (sheets / "combined_biigle_annotations.csv").open(
    newline="", encoding="utf-8-sig"
) as handle:
    data = list(csv.DictReader(handle))

with (sheets / "gl_latlong.csv").open(
    newline="", encoding="utf-8-sig"
) as handle:
    lat_rows = list(csv.DictReader(handle))

with (sheets / "depth_temp_frameid.csv").open(
    newline="", encoding="utf-8-sig"
) as handle:
    meta_rows = list(csv.DictReader(handle))

# -------------------------------------------------------------------------
# Frame/dive universe
# -------------------------------------------------------------------------
frame_to_dive = {}
for row in data:
    frame = clean(row.get("filename"))
    dive = clean(row.get("dive_id"))
    if frame and dive:
        previous = frame_to_dive.get(frame)
        if previous and norm_id(previous) != norm_id(dive):
            report.append((
                "FRAME → dive_id consistency",
                "ERROR",
                f"frame {frame!r} occurs under conflicting dive IDs"
            ))
            fatal = True
        frame_to_dive[frame] = dive

frames = set(frame_to_dive)
dives = {norm_id(x) for x in frame_to_dive.values() if norm_id(x)}

# -------------------------------------------------------------------------
# Latitude / longitude join + duplicate conflicts
# -------------------------------------------------------------------------
coord_lookup = {}
coord_conflicts = []

for line_number, row in enumerate(lat_rows, start=2):
    key = norm_id(row.get("FullID"))
    if not key:
        continue
    value = (clean(row.get("Lat")), clean(row.get("Long")))
    if key in coord_lookup and coord_lookup[key] != value:
        coord_conflicts.append((line_number, key, coord_lookup[key], value))
    else:
        coord_lookup[key] = value

if coord_conflicts:
    report.append((
        "gl_latlong duplicate keys",
        "ERROR",
        f"{len(coord_conflicts)} conflicting duplicate FullID value(s)"
    ))
    fatal = True
else:
    report.append(("gl_latlong duplicate keys", "PASS", "no conflicting duplicate FullID values"))

unmatched_dives = sorted(dives - set(coord_lookup))
report.append((
    "JOIN dive_id ↔ FullID",
    "PASS" if not unmatched_dives else "WARNING",
    f"{len(dives) - len(unmatched_dives)}/{len(dives)} unique dives matched"
))
(out / "00_unmatched_dive_ids.txt").write_text(
    "\n".join(unmatched_dives), encoding="utf-8"
)

# -------------------------------------------------------------------------
# Frame metadata join + duplicate conflicts
# -------------------------------------------------------------------------
meta_lookup = {}
meta_conflicts = []

for line_number, row in enumerate(meta_rows, start=2):
    frame = clean(row.get("FrameID"))
    if not frame:
        continue
    key = frame.casefold()
    values = (
        clean(row.get("video_time")),
        clean(row.get("depth_m")),
        clean(row.get("temperature")),
    )
    if key in meta_lookup and meta_lookup[key] != values:
        meta_conflicts.append((line_number, frame, meta_lookup[key], values))
    else:
        meta_lookup[key] = values

if meta_conflicts:
    report.append((
        "depth_temp_frameid duplicate keys",
        "ERROR",
        f"{len(meta_conflicts)} conflicting duplicate FrameID value(s)"
    ))
    fatal = True
else:
    report.append(("depth_temp_frameid duplicate keys", "PASS", "no conflicting duplicate FrameID values"))

matched_frames = {frame for frame in frames if frame.casefold() in meta_lookup}
unmatched_frames = sorted(frames - matched_frames)

coverage_pct = (100 * len(matched_frames) / len(frames)) if frames else 0.0
report.append((
    "JOIN filename ↔ FrameID",
    "PASS" if not unmatched_frames else "WARNING",
    f"{len(matched_frames)}/{len(frames)} unique BIIGLE frames matched "
    f"({coverage_pct:.1f}% coverage)"
))

with (out / "00_unmatched_frame_ids.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow(["filename"])
    writer.writerows([[x] for x in unmatched_frames])

# -------------------------------------------------------------------------
# Environmental metadata coverage by dive
# -------------------------------------------------------------------------
dive_total = Counter()
dive_matched = Counter()
dive_display = {}

for frame, raw_dive in frame_to_dive.items():
    dive = norm_id(raw_dive)
    dive_total[dive] += 1
    dive_display.setdefault(dive, raw_dive)
    if frame in matched_frames:
        dive_matched[dive] += 1

complete_dives = 0
zero_dives = 0
partial_dives = 0

coverage_rows = []
for dive in sorted(dive_total):
    total = dive_total[dive]
    matched = dive_matched[dive]
    unmatched = total - matched
    pct = 100 * matched / total if total else 0

    if matched == total:
        status = "COMPLETE"
        complete_dives += 1
    elif matched == 0:
        status = "NONE"
        zero_dives += 1
    else:
        status = "PARTIAL"
        partial_dives += 1

    coverage_rows.append([
        dive_display[dive],
        matched,
        total,
        unmatched,
        f"{pct:.1f}",
        status,
    ])

with (out / "00_environmental_metadata_coverage_by_dive.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "dive_id", "matched_frames", "total_frames",
        "unmatched_frames", "coverage_percent", "coverage_status"
    ])
    writer.writerows(coverage_rows)

with (out / "00_environmental_metadata_coverage_summary.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow(["metric", "value"])
    writer.writerows([
        ["total_frames", len(frames)],
        ["matched_frames", len(matched_frames)],
        ["unmatched_frames", len(unmatched_frames)],
        ["frame_coverage_percent", f"{coverage_pct:.1f}"],
        ["total_dives", len(dive_total)],
        ["dives_complete_metadata", complete_dives],
        ["dives_no_metadata", zero_dives],
        ["dives_partial_metadata", partial_dives],
    ])

report.append((
    "Environmental metadata by dive",
    "PASS" if not partial_dives and not zero_dives else "WARNING",
    f"{complete_dives}/{len(dive_total)} dives complete; "
    f"{zero_dives} with no frame metadata; {partial_dives} partially matched"
))

# -------------------------------------------------------------------------
# Raw annotation design QC
# Count each unique Point annotation_id once per frame.
# VME-only points are excluded from the random-point design count.
# Mixed-label points remain one unique point.
# -------------------------------------------------------------------------
annotation_classes = defaultdict(set)
annotation_shape = {}

for row_number, row in enumerate(data, start=2):
    frame = clean(row.get("filename"))
    ann = clean(row.get("annotation_id"))
    if not frame or not ann:
        continue
    key = (frame, ann)
    annotation_classes[key].add(top_level(row.get("label_hierarchy")))
    shape = n(row.get("shape_name"))
    if key in annotation_shape and annotation_shape[key] != shape:
        report.append((
            "annotation_id shape consistency",
            "WARNING",
            f"annotation {ann!r} in {frame!r} has conflicting shape_name values"
        ))
    annotation_shape[key] = shape

point_counts = Counter()
class_conflicts = Counter()

for key, classes in annotation_classes.items():
    frame, ann = key
    if annotation_shape.get(key) != "point":
        continue

    clean_classes = {c for c in classes if c}
    if clean_classes == {"VME"}:
        continue

    point_counts[frame] += 1
    if len(clean_classes) > 1:
        class_conflicts[tuple(sorted(clean_classes))] += 1

frame_point_rows = []
for frame in sorted(frames):
    count = point_counts.get(frame, 0)
    frame_point_rows.append([
        frame,
        frame_to_dive.get(frame, ""),
        count,
        count - TARGET_POINTS,
        "PASS" if count == TARGET_POINTS else "REVIEW",
    ])

with (out / "00_point_design_by_frame.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "filename", "dive_id", "unique_random_point_annotations",
        "difference_from_target", "qc_status"
    ])
    writer.writerows(frame_point_rows)

with (out / "00_frames_not_target_point_count.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "filename", "dive_id", "unique_random_point_annotations",
        "difference_from_target"
    ])
    for row in frame_point_rows:
        if row[2] != TARGET_POINTS:
            writer.writerow(row[:4])

point_sampled_rows = [row for row in frame_point_rows if row[2] > 0]
no_random_point_rows = [row for row in frame_point_rows if row[2] == 0]

exact_target = sum(1 for row in point_sampled_rows if row[2] == TARGET_POINTS)
review_target = len(point_sampled_rows) - exact_target

report.append((
    "Point-sampled frame universe",
    "PASS",
    f"{len(point_sampled_rows)}/{len(frame_point_rows)} raw-export frames contain "
    f"at least one non-VME Point annotation; "
    f"{len(no_random_point_rows)} frame(s) are VME-only/no-random-point frames"
))

report.append((
    f"Random-point design (target {TARGET_POINTS})",
    "PASS" if review_target == 0 else "WARNING",
    f"{exact_target}/{len(point_sampled_rows)} point-sampled frames resolve to exactly "
    f"{TARGET_POINTS} unique non-VME Point annotation IDs; "
    f"{review_target} point-sampled frame(s) flagged for review"
))

with (out / "00_frames_without_random_point_annotations.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "filename", "dive_id", "unique_random_point_annotations",
        "difference_from_target"
    ])
    for row in no_random_point_rows:
        writer.writerow(row[:4])

with (out / "00_annotation_top_level_conflicts.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow(["top_level_combination", "unique_annotation_ids"])
    for combo, count in sorted(class_conflicts.items()):
        writer.writerow([" | ".join(combo), count])

# -------------------------------------------------------------------------
# Final report
# -------------------------------------------------------------------------
with (out / "00_input_validation_summary.csv").open(
    "w", newline="", encoding="utf-8"
) as handle:
    writer = csv.writer(handle)
    writer.writerow(["check", "status", "details"])
    writer.writerows(report)

console_lines = [
    "BIIGLE input validation",
    "=======================",
]
for check, status, details in report:
    console_lines.append(f"{status:7}  {check}: {details}")

console_lines += [
    "",
    "Environmental metadata universe",
    "-------------------------------",
    f"Total BIIGLE frames:            {len(frames):,}",
    f"Frames with depth/temp metadata:{len(matched_frames):9,}",
    f"Frames without metadata:        {len(unmatched_frames):9,}",
    f"Frame metadata coverage:        {coverage_pct:8.1f}%",
    f"Dives with complete metadata:   {complete_dives:9,}/{len(dive_total)}",
    f"Dives with no frame metadata:   {zero_dives:9,}/{len(dive_total)}",
    f"Partially matched dives:        {partial_dives:9,}/{len(dive_total)}",
    "",
    "Random-point design",
    "-------------------",
    f"Target unique points/frame:     {TARGET_POINTS:9,}",
    f"Raw-export frames:              {len(frame_point_rows):9,}",
    f"Point-sampled frames:           {len(point_sampled_rows):9,}",
    f"VME-only/no-point frames:       {len(no_random_point_rows):9,}",
    f"Point frames exactly at target: {exact_target:9,}/{len(point_sampled_rows)}",
    f"Point frames flagged for review:{review_target:9,}/{len(point_sampled_rows)}",
]

console = "\n".join(console_lines) + "\n"
(out / "00_validation_console.txt").write_text(console, encoding="utf-8")
print(console, end="")

if fatal:
    raise SystemExit(1)
PY2
