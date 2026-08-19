#!/usr/bin/env bash

set -euo pipefail

# ============================================================================
# Enrich combined BIIGLE annotations from TAXON_LIST.csv
# ============================================================================
#
# This replaces the OLD generated hierarchy columns:
#   biotic_level_1
#   biotic_level_2
#   biotic_level_3
#   abiotic_level_1
#   abiotic_level_2
#   abiotic_level_3
#
# with the TAXON_LIST fields:
#   cpc_codes
#   kingdom
#   phylum
#   class
#   order
#   family
#   taxonomic_resolution
#   morphology
#   common_id_short
#   common_id_mid
#   common_id_full
#   releif
#   substrate
#   type
#   size
#
# Matching rules:
#   Biotic  -> strip leading "Biotic" from label_hierarchy and match the
#              remaining hierarchy to TAXON_LIST Layer_01 ... Layer_08.
#   Abiotic -> parse relief/substrate/type/size from label_hierarchy and match
#              those fields to the substrate rows in TAXON_LIST.
#   Exclude -> all 15 new fields blank.
#   UNSURE  -> all 15 new fields blank.
#   Other/VME -> blank unless an explicit TAXON_LIST match exists (none is
#                guessed).
#
# If a hierarchy is truncated and matches several TAXON_LIST descendants, only
# values identical across ALL matching TAXON_LIST rows are transferred. Fields
# that would require guessing are left blank.
#
# The six old generated columns are removed. If they already occupy a position
# in the CSV, the 15 new fields are inserted at that same position; otherwise
# the new fields are appended to the end.
#
# The script updates the combined CSV in place after creating a backup.
# It also writes two diagnostic CSVs in the same folder:
#   taxon_lookup_unmatched.csv
#   taxon_lookup_partial_or_ambiguous.csv
#
# Usage:
#   ./extract_hierarchy_columns_TAXON_LIST.sh
#
# Or:
#   ./extract_hierarchy_columns_TAXON_LIST.sh \
#       "/path/to/folder" \
#       "combined_biigle_annotations.csv" \
#       "TAXON_LIST.csv"
# ============================================================================

INPUT_DIR="${1:-.}"
CSV_NAME="${2:-combined_biigle_annotations.csv}"
TAXON_NAME="${3:-TAXON_LIST.csv}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 is required but was not found." >&2
    exit 1
fi

python3 - "$INPUT_DIR" "$CSV_NAME" "$TAXON_NAME" <<'PYTHON'
import csv
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_VERSION = "2026-08-09-taxonomy-substrate-taxON-list-v1"

input_dir = Path(sys.argv[1]).expanduser().resolve()
csv_path = input_dir / sys.argv[2]
taxon_path = input_dir / sys.argv[3]

tmp_path = csv_path.with_name(csv_path.name + ".tmp")
backup_path = csv_path.with_name(
    f"{csv_path.stem}.before_taxon_list_enrichment{csv_path.suffix}"
)
unmatched_report = input_dir / "taxon_lookup_unmatched.csv"
partial_report = input_dir / "taxon_lookup_partial_or_ambiguous.csv"

old_generated_columns = [
    "biotic_level_1",
    "biotic_level_2",
    "biotic_level_3",
    "abiotic_level_1",
    "abiotic_level_2",
    "abiotic_level_3",
]

new_columns = [
    "cpc_codes",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "taxonomic_resolution",
    "morphology",
    "common_id_short",
    "common_id_mid",
    "common_id_full",
    "releif",
    "substrate",
    "type",
    "size",
]

layer_columns = [f"Layer_{i:02d}" for i in range(1, 9)]


def detect_dialect(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        return csv.excel


def normalise_header(value):
    return str(value or "").replace("\ufeff", "").strip().casefold()


def normalise_text(value):
    value = str(value or "").strip().casefold()
    value = value.replace("–", "-").replace("—", "-")
    value = re.sub(r"\s+", " ", value)
    return value


def compact_token(value):
    """Normalise terms used in the Abiotic lookup."""
    value = normalise_text(value)
    value = value.replace("&", "and")
    value = re.sub(r"[\s_/()\-]+", "", value)

    aliases = {
        "wentworthscale": "wentworth",
        "shell(dead)": "shelldead",
        "shelldead": "shelldead",
        "bedformsripples": "ripples",
        "bedformripples": "ripples",
        "ripples": "ripples",
        "bacterialmats": "bacterialmats",
        "bacterialmat": "bacterialmats",
        "icescour": "icescour",
        "bedrockreef": "bedrockreef",
    }

    return aliases.get(value, value)


def find_column(fieldnames, wanted):
    wanted = normalise_header(wanted)
    for fieldname in fieldnames:
        if normalise_header(fieldname) == wanted:
            return fieldname
    return None


def split_hierarchy(value):
    return [
        part.strip()
        for part in str(value or "").split(">")
        if part.strip()
    ]


def hierarchy_key(parts):
    return tuple(normalise_text(part) for part in parts if str(part).strip())


def consensus_values(rows):
    """
    Return only TAXON_LIST values that are identical across all candidate rows.
    Conflicting values are left blank rather than guessed.
    """
    result = {}

    for column in new_columns:
        values = [str(row.get(column, "") or "").strip() for row in rows]
        normalised_nonblank = {
            normalise_text(value): value
            for value in values
            if value != ""
        }

        # A value is safe only if all rows contain the same nonblank value.
        if len(normalised_nonblank) == 1 and all(value != "" for value in values):
            result[column] = next(iter(normalised_nonblank.values()))
        # If every row is blank, keep blank.
        elif all(value == "" for value in values):
            result[column] = ""
        else:
            result[column] = ""

    return result


def row_values(row):
    return {column: str(row.get(column, "") or "").strip() for column in new_columns}


def excel_column_name(number):
    result = ""
    while number:
        number, remainder = divmod(number - 1, 26)
        result = chr(65 + remainder) + result
    return result


def parse_abiotic(parts):
    """
    Parse BIIGLE Abiotic hierarchy into TAXON_LIST substrate dimensions.

    Examples:
      Abiotic > Habitat Relief > Flat > Soft > Wentworth scale > Silt
        -> Flat / soft / wentworth / silt

      Abiotic > Habitat Relief > Flat > Hard > Shell (dead)
        -> Flat / hard / shelldead / [blank]
    """
    values = [normalise_text(part) for part in parts[1:]]

    relief = ""
    substrate = ""
    type_value = ""
    size_value = ""

    if "habitat relief" in values:
        idx = values.index("habitat relief")
        if idx + 1 < len(values):
            relief = values[idx + 1]
    elif values:
        relief = values[0]

    for value in values:
        if value in {"soft", "hard"}:
            substrate = value
            break

    if "wentworth scale" in values:
        idx = values.index("wentworth scale")
        type_value = "wentworth"
        if idx + 1 < len(values):
            size_value = values[idx + 1]
    else:
        # Use the last descriptor after substrate as the type.
        # Examples: Bioturbation, Bacterial mats, Ice scour, Shell (dead).
        if substrate and substrate in values:
            idx = values.index(substrate)
            if idx + 1 < len(values):
                type_value = values[-1]

    return (
        compact_token(relief),
        compact_token(substrate),
        compact_token(type_value),
        compact_token(size_value),
    )


def candidate_abiotic_rows(query_key, abiotic_rows):
    """Return rows consistent with all nonblank parsed Abiotic components."""
    relief, substrate, type_value, size_value = query_key
    candidates = []

    for row, row_key in abiotic_rows:
        keep = True
        for query_component, row_component in zip(query_key, row_key):
            if query_component and query_component != row_component:
                keep = False
                break
        if keep:
            candidates.append(row)

    return candidates


if not input_dir.is_dir():
    raise SystemExit(f"Error: folder does not exist: {input_dir}")

if not csv_path.is_file():
    raise SystemExit(f"Error: combined BIIGLE CSV not found: {csv_path}")

if not taxon_path.is_file():
    raise SystemExit(f"Error: TAXON_LIST.csv not found: {taxon_path}")

# ---------------------------------------------------------------------------
# 1. Read and validate TAXON_LIST.csv
# ---------------------------------------------------------------------------

taxon_dialect = detect_dialect(taxon_path)

with taxon_path.open("r", encoding="utf-8-sig", newline="") as handle:
    taxon_reader = csv.DictReader(handle, dialect=taxon_dialect)

    if not taxon_reader.fieldnames:
        raise SystemExit("Error: TAXON_LIST.csv has no header row.")

    required_taxon_columns = layer_columns + new_columns
    missing_taxon_columns = [
        column
        for column in required_taxon_columns
        if find_column(taxon_reader.fieldnames, column) is None
    ]

    if missing_taxon_columns:
        raise SystemExit(
            "Error: TAXON_LIST.csv is missing required columns: "
            + ", ".join(missing_taxon_columns)
        )

    # Map canonical names to actual source headers.
    taxon_header_map = {
        column: find_column(taxon_reader.fieldnames, column)
        for column in required_taxon_columns
    }

    taxon_rows = []
    for source_row in taxon_reader:
        canonical = {
            column: str(source_row.get(actual, "") or "").strip()
            for column, actual in taxon_header_map.items()
        }
        taxon_rows.append(canonical)

# Taxonomic rows use Layer_01 ... Layer_08.
taxonomy_lookup = defaultdict(list)
taxonomy_rows_with_keys = []

# Substrate rows have blank Layer columns but populated relief/substrate fields.
abiotic_lookup = defaultdict(list)
abiotic_rows_with_keys = []

for row in taxon_rows:
    layers = [row[column] for column in layer_columns if row[column].strip()]

    if layers:
        key = hierarchy_key(layers)
        taxonomy_lookup[key].append(row)
        taxonomy_rows_with_keys.append((row, key))
    elif any(str(row.get(column, "") or "").strip() for column in ["releif", "substrate", "type", "size"]):
        key = (
            compact_token(row.get("releif", "")),
            compact_token(row.get("substrate", "")),
            compact_token(row.get("type", "")),
            compact_token(row.get("size", "")),
        )
        abiotic_lookup[key].append(row)
        abiotic_rows_with_keys.append((row, key))

# ---------------------------------------------------------------------------
# 2. Read combined BIIGLE CSV and establish replacement column position.
# ---------------------------------------------------------------------------

source_dialect = detect_dialect(csv_path)

with csv_path.open("r", encoding="utf-8-sig", newline="") as input_handle:
    reader = csv.DictReader(input_handle, dialect=source_dialect)

    if not reader.fieldnames:
        raise SystemExit("Error: combined BIIGLE CSV has no header row.")

    original_columns = list(reader.fieldnames)
    hierarchy_column = find_column(original_columns, "label_hierarchy")

    if hierarchy_column is None:
        raise SystemExit("Error: could not find label_hierarchy in combined BIIGLE CSV.")

    generated_names = {
        normalise_header(column)
        for column in old_generated_columns + new_columns
    }

    generated_positions = [
        index
        for index, column in enumerate(original_columns)
        if normalise_header(column) in generated_names
    ]

    insertion_original_index = (
        min(generated_positions)
        if generated_positions
        else len(original_columns)
    )

    # Count how many non-generated columns appeared before the insertion point.
    insertion_index = sum(
        1
        for index, column in enumerate(original_columns)
        if index < insertion_original_index
        and normalise_header(column) not in generated_names
    )

    source_columns = [
        column
        for column in original_columns
        if normalise_header(column) not in generated_names
    ]

    output_columns = (
        source_columns[:insertion_index]
        + new_columns
        + source_columns[insertion_index:]
    )

    source_rows = list(reader)

# ---------------------------------------------------------------------------
# 3. Match each BIIGLE label hierarchy to TAXON_LIST.
# ---------------------------------------------------------------------------

counts = Counter()
unmatched_hierarchies = Counter()
partial_hierarchies = Counter()
partial_details = {}

try:
    with tmp_path.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=output_columns)
        writer.writeheader()

        for line_number, row in enumerate(source_rows, start=2):
            output_row = {
                column: row.get(column, "")
                for column in source_columns
            }

            # Always initialise all 15 replacement fields blank.
            for column in new_columns:
                output_row[column] = ""

            hierarchy_text = str(row.get(hierarchy_column, "") or "").strip()
            parts = split_hierarchy(hierarchy_text)
            first_entry = normalise_text(parts[0]) if parts else ""

            # Exclude and UNSURE are never enriched.
            if first_entry == "exclude":
                counts["exclude_blank"] += 1
                writer.writerow(output_row)
                continue

            if first_entry == "unsure":
                counts["unsure_blank"] += 1
                writer.writerow(output_row)
                continue

            matched_values = None
            match_type = None
            candidate_count = 0

            if first_entry == "biotic":
                query_parts = parts[1:]
                key = hierarchy_key(query_parts)

                exact_rows = taxonomy_lookup.get(key, [])

                if exact_rows:
                    candidate_count = len(exact_rows)
                    if len(exact_rows) == 1:
                        matched_values = row_values(exact_rows[0])
                        match_type = "biotic_exact"
                    else:
                        matched_values = consensus_values(exact_rows)
                        match_type = "biotic_exact_ambiguous_consensus"
                elif key:
                    # Truncated hierarchy: use descendants only if supported.
                    prefix_rows = [
                        taxon_row
                        for taxon_row, taxon_key in taxonomy_rows_with_keys
                        if len(taxon_key) >= len(key)
                        and taxon_key[: len(key)] == key
                    ]

                    if prefix_rows:
                        candidate_count = len(prefix_rows)
                        matched_values = consensus_values(prefix_rows)
                        match_type = "biotic_prefix_consensus"

            elif first_entry == "abiotic":
                query_key = parse_abiotic(parts)
                exact_rows = abiotic_lookup.get(query_key, [])

                if exact_rows:
                    candidate_count = len(exact_rows)
                    if len(exact_rows) == 1:
                        matched_values = row_values(exact_rows[0])
                        match_type = "abiotic_exact"
                    else:
                        matched_values = consensus_values(exact_rows)
                        match_type = "abiotic_exact_ambiguous_consensus"
                else:
                    # Incomplete Abiotic hierarchy: preserve only consensus
                    # fields among compatible TAXON_LIST rows.
                    candidate_rows = candidate_abiotic_rows(
                        query_key,
                        abiotic_rows_with_keys,
                    )

                    if candidate_rows:
                        candidate_count = len(candidate_rows)
                        matched_values = consensus_values(candidate_rows)
                        match_type = "abiotic_partial_consensus"

            if matched_values is not None:
                output_row.update(matched_values)
                counts[match_type] += 1

                if "consensus" in match_type:
                    partial_hierarchies[hierarchy_text] += 1
                    partial_details[hierarchy_text] = (
                        match_type,
                        candidate_count,
                    )
            else:
                counts["unmatched"] += 1
                unmatched_hierarchies[hierarchy_text] += 1

            writer.writerow(output_row)

    # -----------------------------------------------------------------------
    # 4. Verify output before replacing original.
    # -----------------------------------------------------------------------

    check_dialect = detect_dialect(tmp_path)
    rows_checked = 0
    exclude_unsure_violations = []

    with tmp_path.open("r", encoding="utf-8-sig", newline="") as check_handle:
        checker = csv.DictReader(check_handle, dialect=check_dialect)
        check_hierarchy = find_column(checker.fieldnames or [], "label_hierarchy")

        missing_output_columns = [
            column
            for column in new_columns
            if find_column(checker.fieldnames or [], column) is None
        ]

        if missing_output_columns:
            raise ValueError(
                "Output verification failed; missing new columns: "
                + ", ".join(missing_output_columns)
            )

        old_columns_remaining = [
            column
            for column in old_generated_columns
            if find_column(checker.fieldnames or [], column) is not None
        ]

        if old_columns_remaining:
            raise ValueError(
                "Output verification failed; old generated columns remain: "
                + ", ".join(old_columns_remaining)
            )

        for line_number, row in enumerate(checker, start=2):
            rows_checked += 1
            parts = split_hierarchy(row.get(check_hierarchy, ""))
            first_entry = normalise_text(parts[0]) if parts else ""

            if first_entry in {"exclude", "unsure"}:
                non_blank = {
                    column: str(row.get(column, "") or "").strip()
                    for column in new_columns
                    if str(row.get(column, "") or "").strip()
                }
                if non_blank:
                    exclude_unsure_violations.append(
                        (line_number, hierarchy_text, non_blank)
                    )

    if rows_checked != len(source_rows):
        raise ValueError(
            f"Row-count verification failed: input={len(source_rows):,}, "
            f"output={rows_checked:,}."
        )

    if exclude_unsure_violations:
        raise ValueError(
            "Verification failed: Exclude/UNSURE rows contain enrichment data."
        )

    # Diagnostic reports.
    with unmatched_report.open("w", encoding="utf-8", newline="") as handle:
        report_writer = csv.writer(handle)
        report_writer.writerow(["label_hierarchy", "row_count"])
        for hierarchy, count in unmatched_hierarchies.most_common():
            report_writer.writerow([hierarchy, count])

    with partial_report.open("w", encoding="utf-8", newline="") as handle:
        report_writer = csv.writer(handle)
        report_writer.writerow(
            ["label_hierarchy", "row_count", "match_type", "candidate_taxon_rows"]
        )
        for hierarchy, count in partial_hierarchies.most_common():
            match_type, candidate_count = partial_details[hierarchy]
            report_writer.writerow(
                [hierarchy, count, match_type, candidate_count]
            )

    shutil.copy2(csv_path, backup_path)
    os.replace(tmp_path, csv_path)

except Exception as error:
    if tmp_path.exists():
        tmp_path.unlink()
    raise SystemExit(f"Error: {error}")

# ---------------------------------------------------------------------------
# 5. Completion summary.
# ---------------------------------------------------------------------------

first_new_position = output_columns.index(new_columns[0]) + 1
last_new_position = output_columns.index(new_columns[-1]) + 1

print(f"\nScript version: {SCRIPT_VERSION}")
print("TAXON_LIST enrichment complete")
print(f"Rows processed:                    {len(source_rows):,}")
print(f"Biotic exact matches:              {counts['biotic_exact']:,}")
print(f"Biotic prefix/consensus matches:   {counts['biotic_prefix_consensus']:,}")
print(f"Biotic ambiguous exact consensus:  {counts['biotic_exact_ambiguous_consensus']:,}")
print(f"Abiotic exact matches:             {counts['abiotic_exact']:,}")
print(f"Abiotic partial/consensus matches: {counts['abiotic_partial_consensus']:,}")
print(f"Exclude rows left blank:           {counts['exclude_blank']:,}")
print(f"UNSURE rows left blank:            {counts['unsure_blank']:,}")
print(f"Unmatched rows left blank:         {counts['unmatched']:,}")
print(f"Updated file:                      {csv_path}")
print(f"Backup file:                       {backup_path}")
print(f"Unmatched report:                  {unmatched_report}")
print(f"Partial/ambiguous report:          {partial_report}")
print(
    "New field block:                    "
    f"{excel_column_name(first_new_position)}:{excel_column_name(last_new_position)}"
)
print("\nNew columns:")
for column in new_columns:
    position = output_columns.index(column) + 1
    print(f"  {excel_column_name(position):>3} ({position:>2}): {column}")
PYTHON
