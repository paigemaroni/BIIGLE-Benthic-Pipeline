#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   ./combine_biigle_csvs.sh [input_folder] [output_csv]
#
# Examples:
#   ./combine_biigle_csvs.sh
#   ./combine_biigle_csvs.sh ./4899_csv_image_annotation_report
#   ./combine_biigle_csvs.sh ./4899_csv_image_annotation_report combined_annotations.csv

INPUT_DIR="${1:-.}"
OUTPUT_CSV="${2:-combined_biigle_annotations.csv}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is required but was not found." >&2
    exit 1
fi

python3 - "$INPUT_DIR" "$OUTPUT_CSV" <<'PYTHON'
import csv
import io
import os
import re
import sys
import zipfile
from pathlib import Path


input_dir = Path(sys.argv[1]).expanduser().resolve()
output_path = Path(sys.argv[2]).expanduser()

# Put a relative output path in the current working directory.
if not output_path.is_absolute():
    output_path = Path.cwd() / output_path

output_path = output_path.resolve()
temporary_output = output_path.with_name(output_path.name + ".tmp")


if not input_dir.is_dir():
    print(
        f"Error: input folder does not exist: {input_dir}",
        file=sys.stderr
    )
    sys.exit(1)


zip_files = sorted(
    path for path in input_dir.glob("*.zip")
    if path.is_file()
)

if not zip_files:
    print(
        f"Error: no ZIP files were found in: {input_dir}",
        file=sys.stderr
    )
    sys.exit(1)


def detect_csv_dialect(zip_archive, member_name):
    """Detect comma-, tab-, or semicolon-delimited files."""

    with zip_archive.open(member_name, "r") as source:
        sample_bytes = source.read(65536)

    sample_text = sample_bytes.decode("utf-8-sig", errors="replace")

    try:
        return csv.Sniffer().sniff(
            sample_text,
            delimiters=",;\t"
        )
    except csv.Error:
        return csv.excel


def extract_dive_id(csv_filename):
    """
    Extract everything beginning with '2025' from the CSV filename.

    Example:
        32633-2025-09-sa-fr01-018.csv
        becomes
        2025-09-sa-fr01-018
    """

    filename_without_extension = Path(csv_filename).stem

    match = re.search(r"(2025.*)$", filename_without_extension)

    if not match:
        raise ValueError(
            f"Could not extract a dive_id beginning with '2025' "
            f"from CSV filename: {csv_filename}"
        )

    return match.group(1)


source_columns = None
writer = None
output_handle = None

zip_count = 0
csv_count = 0
row_count = 0

try:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    for zip_path in zip_files:
        try:
            with zipfile.ZipFile(zip_path, "r") as archive:

                csv_members = sorted(
                    member
                    for member in archive.namelist()
                    if not member.endswith("/")
                    and member.lower().endswith(".csv")
                )

                if not csv_members:
                    print(
                        f"Skipping {zip_path.name}: no CSV file found"
                    )
                    continue

                zip_count += 1

                for member_name in csv_members:

                    # Keep only the CSV filename, excluding any internal
                    # folders within the ZIP archive.
                    parent_csv = Path(member_name).name

                    # Extract the dive ID from the CSV filename.
                    dive_id = extract_dive_id(parent_csv)

                    dialect = detect_csv_dialect(
                        archive,
                        member_name
                    )

                    with archive.open(member_name, "r") as binary_file:
                        with io.TextIOWrapper(
                            binary_file,
                            encoding="utf-8-sig",
                            newline=""
                        ) as text_file:

                            reader = csv.DictReader(
                                text_file,
                                dialect=dialect
                            )

                            current_columns = reader.fieldnames

                            if not current_columns:
                                print(
                                    f"Skipping empty CSV: "
                                    f"{zip_path.name}/{member_name}"
                                )
                                continue

                            # Remove accidental surrounding whitespace
                            # from header names.
                            current_columns = [
                                column.strip() if column else column
                                for column in current_columns
                            ]

                            reader.fieldnames = current_columns

                            if source_columns is None:
                                source_columns = current_columns

                                reserved_columns = {
                                    "dive_id",
                                    "parent_csv"
                                }

                                existing_reserved = (
                                    reserved_columns
                                    & set(source_columns)
                                )

                                if existing_reserved:
                                    raise ValueError(
                                        "A source CSV already contains "
                                        "one or more reserved columns: "
                                        f"{sorted(existing_reserved)}"
                                    )

                                output_handle = temporary_output.open(
                                    "w",
                                    encoding="utf-8",
                                    newline=""
                                )

                                # dive_id is the first column.
                                # parent_csv is the final column.
                                output_columns = (
                                    ["dive_id"]
                                    + source_columns
                                    + ["parent_csv"]
                                )

                                writer = csv.DictWriter(
                                    output_handle,
                                    fieldnames=output_columns,
                                    extrasaction="ignore"
                                )

                                writer.writeheader()

                            elif set(current_columns) != set(source_columns):
                                missing = sorted(
                                    set(source_columns)
                                    - set(current_columns)
                                )

                                additional = sorted(
                                    set(current_columns)
                                    - set(source_columns)
                                )

                                raise ValueError(
                                    f"Column mismatch in "
                                    f"{zip_path.name}/{member_name}\n"
                                    f"Missing columns: {missing}\n"
                                    f"Additional columns: {additional}"
                                )

                            file_row_count = 0

                            for line_number, row in enumerate(
                                reader,
                                start=2
                            ):
                                # DictReader places surplus values under
                                # a None key when a row is malformed.
                                if None in row and row[None]:
                                    raise ValueError(
                                        f"Malformed row in "
                                        f"{zip_path.name}/{member_name} "
                                        f"at line {line_number}: "
                                        f"more values than column names."
                                    )

                                output_row = {
                                    column: row.get(column, "")
                                    for column in source_columns
                                }

                                # Add traceability columns.
                                output_row["dive_id"] = dive_id
                                output_row["parent_csv"] = parent_csv

                                writer.writerow(output_row)

                                file_row_count += 1
                                row_count += 1

                            csv_count += 1

                            print(
                                f"Processed {zip_path.name} -> "
                                f"{parent_csv} -> "
                                f"dive_id={dive_id}: "
                                f"{file_row_count:,} rows"
                            )

        except zipfile.BadZipFile as error:
            raise ValueError(
                f"Invalid or damaged ZIP file: {zip_path}"
            ) from error

    if output_handle is not None:
        output_handle.close()
        output_handle = None

    if csv_count == 0:
        if temporary_output.exists():
            temporary_output.unlink()

        print(
            "Error: ZIP files were found, but none contained "
            "a readable CSV.",
            file=sys.stderr
        )
        sys.exit(1)

    # Replace the final output only after every source has been
    # processed successfully.
    os.replace(temporary_output, output_path)

except Exception as error:
    if output_handle is not None:
        output_handle.close()

    if temporary_output.exists():
        temporary_output.unlink()

    print(f"\nError: {error}", file=sys.stderr)
    sys.exit(1)


print("\nCombination complete")
print(f"ZIP archives processed: {zip_count:,}")
print(f"CSV files processed:     {csv_count:,}")
print(f"Rows written:            {row_count:,}")
print(f"Output file:             {output_path}")
PYTHON