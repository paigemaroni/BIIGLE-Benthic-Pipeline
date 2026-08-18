#!/usr/bin/env bash

set -euo pipefail

# =============================================================================
# SPLIT combined_biigle_annotations.csv INTO POINT DATA + VME DATA
# =============================================================================
#
# Run this script FROM the project folder containing:
#
#   combined_biigle_annotations.csv
#
# The script itself may live inside:
#
#   ./scripts/
#
# Run:
#
#   chmod +x ./scripts/split_point_and_vme_data.sh
#   ./scripts/split_point_and_vme_data.sh
#
# OUTPUTS
# -------
#
#   point_annotations.csv
#       = every row where top_level is NOT VME
#
#   vme_annotations.csv
#       = every row where top_level IS VME
#
# IMPORTANT
# ---------
#
# - The ENTIRE row is copied.
# - ALL original columns are preserved.
# - combined_biigle_annotations.csv is NOT changed.
# - Matching "VME" is case-insensitive and ignores surrounding whitespace.
#
# =============================================================================

INPUT_CSV="combined_biigle_annotations.csv"
POINT_CSV="point_annotations.csv"
VME_CSV="vme_annotations.csv"

PROJECT_DIR="$(pwd)"

INPUT_PATH="${PROJECT_DIR}/${INPUT_CSV}"
POINT_PATH="${PROJECT_DIR}/${POINT_CSV}"
VME_PATH="${PROJECT_DIR}/${VME_CSV}"

if [[ ! -f "$INPUT_PATH" ]]; then
    echo "ERROR: Cannot find ${INPUT_CSV} in:" >&2
    echo "  ${PROJECT_DIR}" >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required." >&2
    exit 1
fi


python3 - "$INPUT_PATH" "$POINT_PATH" "$VME_PATH" <<'PYTHON'

import csv
import sys
from pathlib import Path


input_path = Path(sys.argv[1]).resolve()
point_path = Path(sys.argv[2]).resolve()
vme_path = Path(sys.argv[3]).resolve()


def normalise(value):
    return (
        str(value or "")
        .replace("\ufeff", "")
        .strip()
        .casefold()
    )


def find_column(headers, wanted):
    wanted_norm = normalise(wanted)

    for header in headers:
        if normalise(header) == wanted_norm:
            return header

    return None


total_rows = 0
point_rows = 0
vme_rows = 0


with input_path.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as input_handle, point_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as point_handle, vme_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as vme_handle:

    reader = csv.DictReader(input_handle)

    if not reader.fieldnames:
        raise SystemExit(
            "ERROR: combined_biigle_annotations.csv has no header row."
        )

    top_level_column = find_column(
        reader.fieldnames,
        "top_level",
    )

    if top_level_column is None:
        raise SystemExit(
            "ERROR: Could not find the 'top_level' column."
        )

    fieldnames = list(reader.fieldnames)

    point_writer = csv.DictWriter(
        point_handle,
        fieldnames=fieldnames,
        extrasaction="ignore",
    )

    vme_writer = csv.DictWriter(
        vme_handle,
        fieldnames=fieldnames,
        extrasaction="ignore",
    )

    # Both outputs get the complete original header.
    point_writer.writeheader()
    vme_writer.writeheader()

    for line_number, row in enumerate(
        reader,
        start=2,
    ):
        if None in row and row[None]:
            raise SystemExit(
                f"ERROR: Malformed row at line {line_number}: "
                "more values than column headers."
            )

        total_rows += 1

        if normalise(
            row.get(
                top_level_column,
                "",
            )
        ) == "vme":

            # Write the ENTIRE VME row.
            vme_writer.writerow(row)
            vme_rows += 1

        else:

            # Write the ENTIRE non-VME row.
            point_writer.writerow(row)
            point_rows += 1


# Safety check.
if point_rows + vme_rows != total_rows:
    raise SystemExit(
        "ERROR: Row-count verification failed."
    )


print("")
print("DATA SPLIT COMPLETE")
print("===================")
print(f"Total rows:          {total_rows:,}")
print(f"Point-data rows:     {point_rows:,}")
print(f"VME-data rows:       {vme_rows:,}")
print("")
print("Created:")
print(f"  {point_path}")
print(f"  {vme_path}")
print("")
print("Original file left unchanged:")
print(f"  {input_path}")
print("")
print("Both output files preserve the COMPLETE original rows and columns.")
print("")

PYTHON
