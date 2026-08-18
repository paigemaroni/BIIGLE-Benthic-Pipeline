#!/usr/bin/env bash

set -euo pipefail

# =============================================================================
# SIMPLE FRAME METADATA JOIN
#
# Run from the folder containing:
#   combined_biigle_annotations.csv
#   depth_temp_frameid.csv
#
# This script:
#   1. matches depth_temp_frameid.csv:FrameID
#      to combined_biigle_annotations.csv:filename
#   2. removes any old video_time/depth_m/temperature columns
#   3. removes ONLY trailing unnamed/blank CSV columns
#   4. appends video_time, depth_m, temperature immediately after the
#      true final named column
#   5. keeps unmatched BIIGLE rows with blank values
#
# Usage:
#   chmod +x ./scripts/append_frame_metadata.sh
#   ./scripts/append_frame_metadata.sh
# =============================================================================

MAIN_CSV="combined_biigle_annotations.csv"
LOOKUP_CSV="depth_temp_frameid.csv"

MAIN_PATH="$(pwd)/${MAIN_CSV}"
LOOKUP_PATH="$(pwd)/${LOOKUP_CSV}"

if [[ ! -f "$MAIN_PATH" ]]; then
    echo "ERROR: Cannot find ${MAIN_CSV} in $(pwd)" >&2
    exit 1
fi

if [[ ! -f "$LOOKUP_PATH" ]]; then
    echo "ERROR: Cannot find ${LOOKUP_CSV} in $(pwd)" >&2
    exit 1
fi

python3 - "$MAIN_PATH" "$LOOKUP_PATH" <<'PYTHON'

import csv
import os
import shutil
import sys
from pathlib import Path

main_path = Path(sys.argv[1])
lookup_path = Path(sys.argv[2])

backup_path = main_path.with_name(
    "combined_biigle_annotations.before_frame_metadata.csv"
)
temp_path = main_path.with_name(
    "combined_biigle_annotations.tmp.csv"
)

NEW_COLUMNS = ["video_time", "depth_m", "temperature"]


def norm(value):
    return str(value or "").replace("\ufeff", "").strip().casefold()


def norm_frame(value):
    value = str(value or "").strip().strip('"').strip("'")
    value = value.replace("\\", "/").split("/")[-1].strip()
    return value.casefold()


# -------------------------------------------------------------------------
# 1. Load lookup table
# -------------------------------------------------------------------------

with lookup_path.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)
    try:
        lookup_header = next(reader)
    except StopIteration:
        raise SystemExit("ERROR: depth_temp_frameid.csv is empty.")

    lookup_header_norm = [norm(x) for x in lookup_header]

    def lookup_index(name):
        name = norm(name)
        if name not in lookup_header_norm:
            raise SystemExit(
                f"ERROR: depth_temp_frameid.csv is missing column: {name}"
            )
        return lookup_header_norm.index(name)

    idx_frame = lookup_index("FrameID")
    idx_video = lookup_index("video_time")
    idx_depth = lookup_index("depth_m")
    idx_temp = lookup_index("temperature")

    lookup = {}

    for line_no, row in enumerate(reader, start=2):
        # Pad short rows safely.
        if len(row) < len(lookup_header):
            row = row + [""] * (len(lookup_header) - len(row))

        key = norm_frame(row[idx_frame])

        if not key:
            continue

        values = (
            row[idx_video].strip(),
            row[idx_depth].strip(),
            row[idx_temp].strip(),
        )

        if key in lookup and lookup[key] != values:
            raise SystemExit(
                f"ERROR: conflicting duplicate FrameID at line {line_no}: {key}"
            )

        lookup[key] = values


if not lookup:
    raise SystemExit("ERROR: No usable FrameID values found.")


# -------------------------------------------------------------------------
# 2. Read the main CSV as raw rows
# -------------------------------------------------------------------------

with main_path.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)

    try:
        header = next(reader)
    except StopIteration:
        raise SystemExit("ERROR: combined_biigle_annotations.csv is empty.")

    rows = list(reader)


# Make every row the same width as the header before editing.
width = len(header)

for i, row in enumerate(rows):
    if len(row) < width:
        rows[i] = row + [""] * (width - len(row))
    elif len(row) > width:
        # Preserve extra data by extending the header with blank names.
        extra = len(row) - width
        header.extend([""] * extra)
        width = len(header)

        for j in range(i):
            if len(rows[j]) < width:
                rows[j].extend([""] * (width - len(rows[j])))


# -------------------------------------------------------------------------
# 3. Remove old generated columns wherever they currently occur
# -------------------------------------------------------------------------

remove_names = {norm(x) for x in NEW_COLUMNS}

keep_indices = [
    i
    for i, name in enumerate(header)
    if norm(name) not in remove_names
]

header = [header[i] for i in keep_indices]
rows = [
    [row[i] if i < len(row) else "" for i in keep_indices]
    for row in rows
]


# -------------------------------------------------------------------------
# 4. Remove ONLY trailing unnamed columns
#
# This is the key fix for columns jumping out to FUQ/FUR/etc.
# Named columns are preserved even if their data cells are blank.
# -------------------------------------------------------------------------

last_named_index = -1

for i, name in enumerate(header):
    if str(name or "").strip() != "":
        last_named_index = i

if last_named_index < 0:
    raise SystemExit("ERROR: No named columns found in main CSV.")

header = header[: last_named_index + 1]
rows = [
    row[: last_named_index + 1]
    for row in rows
]


# -------------------------------------------------------------------------
# 5. Find filename column
# -------------------------------------------------------------------------

header_norm = [norm(x) for x in header]

if "filename" not in header_norm:
    raise SystemExit(
        "ERROR: Could not find filename column in combined_biigle_annotations.csv."
    )

filename_index = header_norm.index("filename")


# -------------------------------------------------------------------------
# 6. Append new columns immediately after the true last named column
# -------------------------------------------------------------------------

output_header = header + NEW_COLUMNS

matched_rows = 0
matched_frames = set()
unmatched_rows = 0

output_rows = []

for row in rows:
    if len(row) < len(header):
        row = row + [""] * (len(header) - len(row))

    key = norm_frame(row[filename_index])

    if key in lookup:
        video_time, depth_m, temperature = lookup[key]
        matched_rows += 1
        matched_frames.add(key)
    else:
        video_time = ""
        depth_m = ""
        temperature = ""
        unmatched_rows += 1

    output_rows.append(
        row + [video_time, depth_m, temperature]
    )


if matched_rows == 0:
    raise SystemExit(
        "ERROR: Zero FrameID/filename matches found. Original file unchanged."
    )


# -------------------------------------------------------------------------
# 7. Write temporary file, back up original, replace original
# -------------------------------------------------------------------------

with temp_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(output_header)
    writer.writerows(output_rows)

shutil.copy2(main_path, backup_path)
os.replace(temp_path, main_path)


# -------------------------------------------------------------------------
# 8. Report exact new column positions
# -------------------------------------------------------------------------

def excel_col(n):
    result = ""
    while n:
        n, rem = divmod(n - 1, 26)
        result = chr(65 + rem) + result
    return result


first_new = len(header) + 1

print("")
print("SUCCESS")
print("=======")
print(f"Lookup FrameIDs loaded:   {len(lookup):,}")
print(f"Unique filenames matched: {len(matched_frames):,}")
print(f"BIIGLE rows updated:      {matched_rows:,}")
print(f"BIIGLE rows left blank:   {unmatched_rows:,}")
print("")
print("New columns placed immediately after the final named column:")
for offset, name in enumerate(NEW_COLUMNS):
    col_num = first_new + offset
    print(f"  {excel_col(col_num)} (column {col_num}): {name}")
print("")
print(f"Updated: {main_path}")
print(f"Backup:  {backup_path}")
print("")

PYTHON
