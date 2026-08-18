#!/usr/bin/env bash

set -euo pipefail

# Add or replace a "top_level" column in combined_biigle_annotations.csv.
#
# The value is taken from everything before the first ">"
# in the label_hierarchy column.
#
# The top_level column is inserted as COLUMN D (the fourth column).
# All columns that previously appeared from column D onward move one place right.
#
# Examples:
#   "Animalia > Porifera > Demospongiae" -> "Animalia"
#   "Porifera"                           -> "Porifera"
#   ""                                   -> ""
#
# The script updates the CSV in place and creates:
#   combined_biigle_annotations.before_top_level.csv
#
# Usage:
#   ./add_top_level_column_D.sh
#
# Or:
#   ./add_top_level_column_D.sh "/path/to/folder"
#
# Optional filename:
#   ./add_top_level_column_D.sh \
#       "/path/to/folder" \
#       "combined_biigle_annotations.csv"

INPUT_DIR="${1:-.}"
CSV_NAME="${2:-combined_biigle_annotations.csv}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is required but was not found." >&2
    exit 1
fi

python3 - "$INPUT_DIR" "$CSV_NAME" <<'PYTHON'

import csv
import os
import shutil
import sys
from pathlib import Path


input_dir = Path(sys.argv[1]).expanduser().resolve()
csv_path = input_dir / sys.argv[2]

backup_path = csv_path.with_name(
    f"{csv_path.stem}.before_top_level{csv_path.suffix}"
)
temporary_path = csv_path.with_name(
    f"{csv_path.name}.tmp"
)


def detect_dialect(path):
    """Detect comma-, tab-, or semicolon-delimited files."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)

    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        return csv.excel


def locate_column(fieldnames, target):
    """Find a column case-insensitively, ignoring surrounding spaces."""
    target = target.strip().casefold()

    for fieldname in fieldnames:
        if (fieldname or "").strip().casefold() == target:
            return fieldname

    return None


if not input_dir.is_dir():
    print(f"Error: folder does not exist: {input_dir}", file=sys.stderr)
    sys.exit(1)

if not csv_path.is_file():
    print(f"Error: CSV file not found: {csv_path}", file=sys.stderr)
    sys.exit(1)


dialect = detect_dialect(csv_path)

row_count = 0
blank_count = 0

try:
    with csv_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as input_handle, temporary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as output_handle:

        reader = csv.DictReader(input_handle, dialect=dialect)

        if not reader.fieldnames:
            raise ValueError("The CSV does not contain a header row.")

        hierarchy_column = locate_column(
            reader.fieldnames,
            "label_hierarchy",
        )

        if hierarchy_column is None:
            raise ValueError(
                "Could not find a 'label_hierarchy' column. "
                "Columns found: "
                + ", ".join(repr(name) for name in reader.fieldnames)
            )

        existing_top_level = locate_column(
            reader.fieldnames,
            "top_level",
        )

        # Remove any existing top_level column so the script can be rerun.
        source_columns = [
            column
            for column in reader.fieldnames
            if column != existing_top_level
        ]

        if len(source_columns) < 3:
            raise ValueError(
                "The CSV has fewer than three existing columns, "
                "so top_level cannot be inserted as column D."
            )

        # Insert top_level at zero-based index 3:
        # A = 0, B = 1, C = 2, D = 3.
        output_columns = (
            source_columns[:3]
            + ["top_level"]
            + source_columns[3:]
        )

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

            hierarchy = str(
                row.get(hierarchy_column, "") or ""
            ).strip()

            if hierarchy:
                # Split only at the first hierarchy separator.
                top_level = hierarchy.split(">", 1)[0].strip()
            else:
                top_level = ""
                blank_count += 1

            output_row = {
                column: row.get(column, "")
                for column in source_columns
            }
            output_row["top_level"] = top_level

            writer.writerow(output_row)
            row_count += 1

    # Preserve the source CSV before replacing it.
    shutil.copy2(csv_path, backup_path)
    os.replace(temporary_path, csv_path)

except Exception as error:
    if temporary_path.exists():
        temporary_path.unlink()

    print(f"Error: {error}", file=sys.stderr)
    sys.exit(1)


print("\nTop-level extraction complete")
print(f"Rows processed:               {row_count:,}")
print(f"Blank label_hierarchy rows:   {blank_count:,}")
print("top_level position:           Column D")
print(f"Updated file:                 {csv_path}")
print(f"Backup file:                  {backup_path}")

PYTHON
