#!/usr/bin/env bash

set -euo pipefail

# Robustly add Lat/Long from gl_latlong.csv to combined_biigle_annotations.csv.
#
# Matching:
#   gl_latlong.csv                  FullID
#   combined_biigle_annotations.csv dive_id
#
# The script normalises both IDs by:
#   - removing folder paths
#   - removing a trailing .csv
#   - extracting the portion beginning with 2025
#   - trimming whitespace
#   - matching case-insensitively
#
# It updates combined_biigle_annotations.csv IN PLACE and first creates:
#   combined_biigle_annotations.backup.csv
#
# Unmatched rows are retained with blank lat and long values.
#
# Usage:
#   ./add_latlong_to_biigle_v2.sh
#
# Or:
#   ./add_latlong_to_biigle_v2.sh "/path/to/folder"
#
# Optional explicit filenames:
#   ./add_latlong_to_biigle_v2.sh \
#       "/path/to/folder" \
#       "combined_biigle_annotations.csv" \
#       "gl_latlong.csv"

INPUT_DIR="${1:-.}"
COMBINED_NAME="${2:-combined_biigle_annotations.csv}"
LOOKUP_NAME="${3:-gl_latlong.csv}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is required but was not found." >&2
    exit 1
fi

python3 - "$INPUT_DIR" "$COMBINED_NAME" "$LOOKUP_NAME" <<'PYTHON'

import csv
import os
import re
import shutil
import sys
from pathlib import Path


input_dir = Path(sys.argv[1]).expanduser().resolve()
combined_path = input_dir / sys.argv[2]
lookup_path = input_dir / sys.argv[3]

backup_path = combined_path.with_name(
    f"{combined_path.stem}.backup{combined_path.suffix}"
)
temporary_path = combined_path.with_name(
    f"{combined_path.name}.tmp"
)


def detect_dialect(path):
    """Detect comma-, tab-, or semicolon-delimited files."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)

    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        return csv.excel


def clean_header(value):
    """Normalise a header for case-insensitive matching."""
    return (value or "").replace("\ufeff", "").strip().casefold()


def locate_column(fieldnames, accepted_names):
    """Find a column using case-insensitive, whitespace-tolerant matching."""
    accepted = {clean_header(name) for name in accepted_names}

    for actual_name in fieldnames:
        if clean_header(actual_name) in accepted:
            return actual_name

    return None


def normalise_id(value):
    """
    Convert FullID or dive_id into a common matching key.

    Examples:
      32633-2025-09-sa-fr01-018.csv -> 2025-09-sa-fr01-018
      folder/2025-09-sa-fr01-018.csv -> 2025-09-sa-fr01-018
      2025-09-SA-FR01-018            -> 2025-09-sa-fr01-018
    """
    value = str(value or "").strip().strip('"').strip("'")

    # Treat both slash types as path separators.
    value = value.replace("\\", "/").split("/")[-1].strip()

    # Remove a trailing CSV extension.
    value = re.sub(r"\.csv$", "", value, flags=re.IGNORECASE).strip()

    # Extract everything beginning with 2025.
    match = re.search(r"2025.*$", value, flags=re.IGNORECASE)
    if match:
        value = match.group(0)

    # Normalise common dash characters and whitespace.
    value = (
        value.replace("–", "-")
             .replace("—", "-")
             .replace("_", "-")
    )
    value = re.sub(r"\s+", "", value)

    return value.casefold()


if not input_dir.is_dir():
    print(f"Error: folder does not exist: {input_dir}", file=sys.stderr)
    sys.exit(1)

if not combined_path.is_file():
    print(
        f"Error: combined annotation file not found: {combined_path}",
        file=sys.stderr,
    )
    sys.exit(1)

if not lookup_path.is_file():
    print(
        f"Error: coordinate file not found: {lookup_path}",
        file=sys.stderr,
    )
    sys.exit(1)


# -------------------------------------------------------------------------
# Read gl_latlong.csv and create an ID -> coordinates lookup.
# -------------------------------------------------------------------------

lookup_dialect = detect_dialect(lookup_path)
coordinate_lookup = {}
lookup_display_ids = {}
duplicate_identical = 0
conflicts = []

with lookup_path.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:
    reader = csv.DictReader(handle, dialect=lookup_dialect)

    if not reader.fieldnames:
        print("Error: gl_latlong.csv has no header.", file=sys.stderr)
        sys.exit(1)

    full_id_column = locate_column(
        reader.fieldnames,
        {"FullID", "Full ID", "full_id"},
    )
    lat_column = locate_column(
        reader.fieldnames,
        {"Lat", "Latitude"},
    )
    long_column = locate_column(
        reader.fieldnames,
        {"Long", "Longitude", "Lon", "Lng"},
    )

    missing = []
    if full_id_column is None:
        missing.append("FullID")
    if lat_column is None:
        missing.append("Lat")
    if long_column is None:
        missing.append("Long")

    if missing:
        print(
            "Error: gl_latlong.csv is missing required columns: "
            + ", ".join(missing),
            file=sys.stderr,
        )
        print(
            "Columns actually found: "
            + ", ".join(repr(name) for name in reader.fieldnames),
            file=sys.stderr,
        )
        sys.exit(1)

    for line_number, row in enumerate(reader, start=2):
        raw_full_id = row.get(full_id_column, "")
        match_id = normalise_id(raw_full_id)
        lat = str(row.get(lat_column, "") or "").strip()
        long = str(row.get(long_column, "") or "").strip()

        if not match_id:
            continue

        new_coordinates = (lat, long)

        if match_id in coordinate_lookup:
            if coordinate_lookup[match_id] == new_coordinates:
                duplicate_identical += 1
            else:
                conflicts.append(
                    (
                        raw_full_id,
                        coordinate_lookup[match_id],
                        new_coordinates,
                        line_number,
                    )
                )
        else:
            coordinate_lookup[match_id] = new_coordinates
            lookup_display_ids[match_id] = str(raw_full_id)


if conflicts:
    print(
        "Error: duplicate FullID values have conflicting coordinates.",
        file=sys.stderr,
    )
    for raw_id, first_coords, second_coords, line_number in conflicts[:20]:
        print(
            f"  {raw_id!r}: {first_coords} versus "
            f"{second_coords} at line {line_number}",
            file=sys.stderr,
        )
    sys.exit(1)

if not coordinate_lookup:
    print(
        "Error: no usable FullID coordinate records were read from "
        "gl_latlong.csv.",
        file=sys.stderr,
    )
    sys.exit(1)


# -------------------------------------------------------------------------
# Join coordinates to the combined annotation file.
# -------------------------------------------------------------------------

combined_dialect = detect_dialect(combined_path)

rows_total = 0
rows_matched = 0
rows_unmatched = 0
unique_matched = set()
unique_unmatched = set()
sample_combined_ids = []
sample_lookup_ids = list(lookup_display_ids.values())[:10]

try:
    with combined_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as input_handle, temporary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as output_handle:

        reader = csv.DictReader(input_handle, dialect=combined_dialect)

        if not reader.fieldnames:
            raise ValueError(
                "combined_biigle_annotations.csv has no header."
            )

        dive_id_column = locate_column(reader.fieldnames, {"dive_id", "Dive ID"})

        if dive_id_column is None:
            raise ValueError(
                "Could not find the dive_id column. Columns found: "
                + ", ".join(repr(name) for name in reader.fieldnames)
            )

        # Remove existing lat/long columns so the script can be rerun safely.
        existing_lat_column = locate_column(
            reader.fieldnames,
            {"lat", "latitude"},
        )
        existing_long_column = locate_column(
            reader.fieldnames,
            {"long", "longitude", "lon", "lng"},
        )

        source_columns = [
            column
            for column in reader.fieldnames
            if column not in {existing_lat_column, existing_long_column}
        ]

        output_columns = source_columns + ["lat", "long"]

        writer = csv.DictWriter(
            output_handle,
            fieldnames=output_columns,
            extrasaction="ignore",
        )
        writer.writeheader()

        for line_number, row in enumerate(reader, start=2):
            if None in row and row[None]:
                raise ValueError(
                    f"Malformed row at line {line_number}: "
                    "more values than headers."
                )

            rows_total += 1

            raw_dive_id = row.get(dive_id_column, "")
            match_id = normalise_id(raw_dive_id)

            if len(sample_combined_ids) < 10 and raw_dive_id:
                sample_combined_ids.append(str(raw_dive_id))

            if match_id and match_id in coordinate_lookup:
                lat, long = coordinate_lookup[match_id]
                rows_matched += 1
                unique_matched.add(match_id)
            else:
                # Keep unmatched rows, but leave coordinates blank.
                lat, long = "", ""
                rows_unmatched += 1
                if match_id:
                    unique_unmatched.add(match_id)

            output_row = {
                column: row.get(column, "")
                for column in source_columns
            }
            output_row["lat"] = lat
            output_row["long"] = long
            writer.writerow(output_row)

    # Make a backup of the current combined file before replacing it.
    shutil.copy2(combined_path, backup_path)
    os.replace(temporary_path, combined_path)

except Exception as error:
    if temporary_path.exists():
        temporary_path.unlink()

    print(f"Error: {error}", file=sys.stderr)
    sys.exit(1)


# -------------------------------------------------------------------------
# Results and diagnostics.
# -------------------------------------------------------------------------

print("\nLatitude/longitude join complete")
print(f"Lookup IDs loaded:          {len(coordinate_lookup):,}")
print(f"Identical duplicate IDs:    {duplicate_identical:,}")
print(f"Annotation rows processed:  {rows_total:,}")
print(f"Annotation rows matched:    {rows_matched:,}")
print(f"Annotation rows unmatched:  {rows_unmatched:,}")
print(f"Unique dive IDs matched:    {len(unique_matched):,}")
print(f"Unique dive IDs unmatched:  {len(unique_unmatched):,}")
print(f"Updated file:               {combined_path}")
print(f"Backup file:                {backup_path}")

if rows_matched == 0:
    print(
        "\nWARNING: zero rows matched. Here are example values "
        "to help diagnose the mismatch:"
    )

    print("\nExample dive_id values from combined annotations:")
    for value in sample_combined_ids:
        print(f"  raw={value!r}  normalised={normalise_id(value)!r}")

    print("\nExample FullID values from gl_latlong.csv:")
    for value in sample_lookup_ids:
        print(f"  raw={value!r}  normalised={normalise_id(value)!r}")

elif unique_unmatched:
    print("\nUnmatched dive IDs were retained with blank lat/long:")
    for match_id in sorted(unique_unmatched)[:100]:
        print(f"  - {match_id}")

    if len(unique_unmatched) > 100:
        print(
            f"  ...and {len(unique_unmatched) - 100:,} more."
        )

PYTHON
