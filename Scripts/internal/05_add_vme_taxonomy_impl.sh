#!/usr/bin/env bash

set -euo pipefail

# =============================================================================
# ADD VME TAXONOMY TO combined_biigle_annotations.csv
# =============================================================================
#
# Run this script FROM the project folder containing:
#
#   combined_biigle_annotations.csv
#   vme_taxon_list.csv
#
# The shell script itself may live inside ./scripts/
#
# Example:
#
#   chmod +x ./scripts/add_vme_taxonomy.sh
#   ./scripts/add_vme_taxonomy.sh
#
# WHAT IT DOES
# ------------
#
# Only rows where:
#
#   top_level == VME
#
# are modified.
#
# It matches the BIIGLE hierarchy:
#
#   VME > Layer_01 > Layer_02 > Layer_03 > Layer_04 > Layer_05
#
# against:
#
#   vme_taxon_list.csv
#
# and fills these columns:
#
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
#
# Existing Biotic, Abiotic, Exclude and UNSURE rows are NOT changed.
#
# MATCHING
# --------
#
# 1. Exact hierarchy match is attempted first.
#
# 2. If exact matching fails, a conservative ordered-layer match is used.
#    This allows BIIGLE to contain an additional intermediate taxonomic level
#    that is omitted from the lookup table.
#
#    Example:
#
#    BIIGLE:
#      VME > Tunicates (...) > Tunicata > Ascidiacea >
#      Dense solitary ascidian aggregation
#
#    Lookup:
#      Tunicates (...) > Tunicata > Dense solitary ascidian aggregation
#
#    The lookup layers must still occur in the same order and the best match
#    must be unique. Ambiguous matches are NOT guessed.
#
# BACKUP
# ------
#
# Before replacing the main CSV, the script creates:
#
#   combined_biigle_annotations.before_vme_taxonomy.csv
#
# REPORTS
# -------
#
#   vme_taxonomy_unmatched.csv
#   vme_taxonomy_ambiguous.csv
#
# =============================================================================


MAIN_CSV="combined_biigle_annotations.csv"
LOOKUP_CSV="vme_taxon_list.csv"

PROJECT_DIR="$(pwd)"
MAIN_PATH="${PROJECT_DIR}/${MAIN_CSV}"
LOOKUP_PATH="${PROJECT_DIR}/${LOOKUP_CSV}"

if [[ ! -f "$MAIN_PATH" ]]; then
    echo "ERROR: Cannot find ${MAIN_CSV} in:" >&2
    echo "  ${PROJECT_DIR}" >&2
    echo "" >&2
    echo "Run this script from the project folder containing the CSV." >&2
    exit 1
fi

if [[ ! -f "$LOOKUP_PATH" ]]; then
    echo "ERROR: Cannot find ${LOOKUP_CSV} in:" >&2
    echo "  ${PROJECT_DIR}" >&2
    echo "" >&2
    echo "Run this script from the project folder containing the CSV." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required." >&2
    exit 1
fi


python3 - "$MAIN_PATH" "$LOOKUP_PATH" <<'PYTHON'

import csv
import os
import shutil
import sys
from pathlib import Path


main_path = Path(sys.argv[1]).resolve()
lookup_path = Path(sys.argv[2]).resolve()

backup_path = main_path.with_name(
    "combined_biigle_annotations.before_vme_taxonomy.csv"
)

temp_path = main_path.with_name(
    "combined_biigle_annotations.vme_taxonomy_tmp.csv"
)

unmatched_path = main_path.with_name(
    "vme_taxonomy_unmatched.csv"
)

ambiguous_path = main_path.with_name(
    "vme_taxonomy_ambiguous.csv"
)


TAXON_COLUMNS = [
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
]

LAYER_COLUMNS = [
    "Layer_01",
    "Layer_02",
    "Layer_03",
    "Layer_04",
    "Layer_05",
]


def normalise(value):
    """
    Normalise text for matching while preserving the original values
    written back to the dataset.
    """
    return " ".join(
        str(value or "")
        .replace("\ufeff", "")
        .strip()
        .casefold()
        .split()
    )


def find_column(headers, wanted):
    wanted_norm = normalise(wanted)

    for header in headers:
        if normalise(header) == wanted_norm:
            return header

    return None


def split_hierarchy(value):
    return [
        part.strip()
        for part in str(value or "").split(">")
        if part.strip()
    ]


def is_ordered_subsequence(lookup_layers, hierarchy_layers):
    """
    True when every nonblank lookup layer occurs in the hierarchy
    in the same order, allowing extra hierarchy levels between them.
    """
    position = 0

    for layer in lookup_layers:
        found = False

        while position < len(hierarchy_layers):
            if hierarchy_layers[position] == layer:
                found = True
                position += 1
                break

            position += 1

        if not found:
            return False

    return True


def candidate_score(lookup_layers, hierarchy_layers):
    """
    Score a conservative fallback candidate.

    Priority:
      1. final lookup layer equals final BIIGLE hierarchy layer;
      2. more lookup layers = more specific;
      3. more adjacent transitions retained in BIIGLE.
    """
    if not lookup_layers:
        return (-1, -1, -1)

    final_equal = int(
        bool(hierarchy_layers)
        and lookup_layers[-1] == hierarchy_layers[-1]
    )

    positions = []
    start = 0

    for layer in lookup_layers:
        found_position = None

        for index in range(start, len(hierarchy_layers)):
            if hierarchy_layers[index] == layer:
                found_position = index
                break

        if found_position is None:
            return (-1, -1, -1)

        positions.append(found_position)
        start = found_position + 1

    adjacent_transitions = sum(
        1
        for first, second in zip(
            positions,
            positions[1:],
        )
        if second == first + 1
    )

    return (
        final_equal,
        len(lookup_layers),
        adjacent_transitions,
    )


# =============================================================================
# 1. READ VME TAXON LOOKUP
# =============================================================================

lookup_rows = []
exact_lookup = {}
duplicate_keys = []

with lookup_path.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:

    reader = csv.DictReader(handle)

    if not reader.fieldnames:
        raise SystemExit(
            "ERROR: vme_taxon_list.csv has no header row."
        )

    layer_headers = {}

    for required in LAYER_COLUMNS:
        actual = find_column(
            reader.fieldnames,
            required,
        )

        if actual is None:
            raise SystemExit(
                f"ERROR: vme_taxon_list.csv is missing {required}."
            )

        layer_headers[required] = actual

    taxon_headers = {}

    for required in TAXON_COLUMNS:
        actual = find_column(
            reader.fieldnames,
            required,
        )

        if actual is None:
            raise SystemExit(
                f"ERROR: vme_taxon_list.csv is missing {required}."
            )

        taxon_headers[required] = actual

    for line_number, row in enumerate(
        reader,
        start=2,
    ):
        layers_original = [
            str(
                row.get(
                    layer_headers[layer],
                    "",
                )
                or ""
            ).strip()
            for layer in LAYER_COLUMNS
        ]

        layers_normalised = [
            normalise(value)
            for value in layers_original
        ]

        nonblank_layers = [
            value
            for value in layers_normalised
            if value
        ]

        if not nonblank_layers:
            continue

        exact_key = tuple(
            layers_normalised
        )

        taxon_values = {
            target: str(
                row.get(
                    taxon_headers[target],
                    "",
                )
                or ""
            ).strip()
            for target in TAXON_COLUMNS
        }

        lookup_record = {
            "line_number": line_number,
            "layers_original": layers_original,
            "layers_normalised": layers_normalised,
            "nonblank_layers": nonblank_layers,
            "taxon_values": taxon_values,
        }

        if exact_key in exact_lookup:
            duplicate_keys.append(
                (
                    line_number,
                    layers_original,
                )
            )
        else:
            exact_lookup[exact_key] = (
                lookup_record
            )

        lookup_rows.append(
            lookup_record
        )


if duplicate_keys:
    raise SystemExit(
        "ERROR: Duplicate Layer_01..Layer_05 combinations exist "
        "in vme_taxon_list.csv."
    )

if not lookup_rows:
    raise SystemExit(
        "ERROR: No usable VME lookup rows were found."
    )


# =============================================================================
# 2. MATCH FUNCTION
# =============================================================================

def match_vme_hierarchy(label_hierarchy):
    hierarchy_parts = split_hierarchy(
        label_hierarchy
    )

    if (
        hierarchy_parts
        and normalise(hierarchy_parts[0]) == "vme"
    ):
        hierarchy_parts = hierarchy_parts[1:]

    hierarchy_normalised = [
        normalise(value)
        for value in hierarchy_parts
    ]

    # -------------------------------------------------------------
    # Exact match first.
    # Pad/truncate BIIGLE hierarchy to five lookup layers.
    # -------------------------------------------------------------

    exact_key = tuple(
        hierarchy_normalised[index]
        if index < len(hierarchy_normalised)
        else ""
        for index in range(5)
    )

    if exact_key in exact_lookup:
        return (
            "exact",
            exact_lookup[exact_key],
            [],
        )

    # -------------------------------------------------------------
    # Conservative fallback.
    # -------------------------------------------------------------

    candidates = []

    for lookup_record in lookup_rows:
        lookup_layers = lookup_record[
            "nonblank_layers"
        ]

        if not hierarchy_normalised:
            continue

        # Layer_01 must agree.
        if (
            not lookup_layers
            or lookup_layers[0]
            != hierarchy_normalised[0]
        ):
            continue

        if not is_ordered_subsequence(
            lookup_layers,
            hierarchy_normalised,
        ):
            continue

        score = candidate_score(
            lookup_layers,
            hierarchy_normalised,
        )

        candidates.append(
            (
                score,
                lookup_record,
            )
        )

    if not candidates:
        return (
            "unmatched",
            None,
            [],
        )

    best_score = max(
        score
        for score, record in candidates
    )

    best_candidates = [
        record
        for score, record in candidates
        if score == best_score
    ]

    if len(best_candidates) == 1:
        return (
            "fallback",
            best_candidates[0],
            [],
        )

    return (
        "ambiguous",
        None,
        best_candidates,
    )


# =============================================================================
# 3. READ AND UPDATE combined_biigle_annotations.csv
# =============================================================================

exact_matches = 0
fallback_matches = 0
unmatched_rows = []
ambiguous_rows = []
vme_rows = 0
non_vme_rows = 0

try:
    with main_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as input_handle, temp_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as output_handle:

        reader = csv.DictReader(
            input_handle
        )

        if not reader.fieldnames:
            raise ValueError(
                "combined_biigle_annotations.csv has no header row."
            )

        top_level_column = find_column(
            reader.fieldnames,
            "top_level",
        )

        hierarchy_column = find_column(
            reader.fieldnames,
            "label_hierarchy",
        )

        if top_level_column is None:
            raise ValueError(
                "Could not find top_level column."
            )

        if hierarchy_column is None:
            raise ValueError(
                "Could not find label_hierarchy column."
            )

        # ---------------------------------------------------------
        # Preserve every existing column in its current order.
        # If any of the 11 target columns are missing, append only
        # the missing ones to the end.
        # ---------------------------------------------------------

        output_columns = list(
            reader.fieldnames
        )

        target_column_map = {}

        for target in TAXON_COLUMNS:
            existing = find_column(
                output_columns,
                target,
            )

            if existing is None:
                output_columns.append(
                    target
                )
                target_column_map[target] = (
                    target
                )
            else:
                target_column_map[target] = (
                    existing
                )

        writer = csv.DictWriter(
            output_handle,
            fieldnames=output_columns,
            extrasaction="ignore",
        )

        writer.writeheader()

        for line_number, row in enumerate(
            reader,
            start=2,
        ):
            # Ensure any newly appended target fields exist.
            output_row = {
                column: row.get(
                    column,
                    "",
                )
                for column in output_columns
            }

            if (
                normalise(
                    row.get(
                        top_level_column,
                        "",
                    )
                )
                != "vme"
            ):
                non_vme_rows += 1
                writer.writerow(
                    output_row
                )
                continue

            vme_rows += 1

            hierarchy_value = row.get(
                hierarchy_column,
                "",
            )

            match_type, lookup_record, candidates = (
                match_vme_hierarchy(
                    hierarchy_value
                )
            )

            if match_type in {
                "exact",
                "fallback",
            }:
                # Overwrite only the 11 VME taxonomic fields.
                # All other columns remain untouched.
                for target in TAXON_COLUMNS:
                    output_column = (
                        target_column_map[target]
                    )

                    output_row[
                        output_column
                    ] = lookup_record[
                        "taxon_values"
                    ][target]

                if match_type == "exact":
                    exact_matches += 1
                else:
                    fallback_matches += 1

            elif match_type == "unmatched":
                unmatched_rows.append(
                    {
                        "csv_line": line_number,
                        "label_hierarchy": hierarchy_value,
                    }
                )

            else:
                ambiguous_rows.append(
                    {
                        "csv_line": line_number,
                        "label_hierarchy": hierarchy_value,
                        "candidate_lookup_lines": " | ".join(
                            str(
                                candidate[
                                    "line_number"
                                ]
                            )
                            for candidate
                            in candidates
                        ),
                        "candidate_layers": " || ".join(
                            " > ".join(
                                value
                                for value in candidate[
                                    "layers_original"
                                ]
                                if value
                            )
                            for candidate
                            in candidates
                        ),
                    }
                )

            writer.writerow(
                output_row
            )

    # -------------------------------------------------------------
    # Safety verification before replacing the main file.
    # -------------------------------------------------------------

    with temp_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as check_handle:

        checker = csv.DictReader(
            check_handle
        )

        if not checker.fieldnames:
            raise ValueError(
                "Output verification failed: no header row."
            )

        for target in TAXON_COLUMNS:
            if find_column(
                checker.fieldnames,
                target,
            ) is None:
                raise ValueError(
                    "Output verification failed: "
                    f"missing {target}."
                )

    # -------------------------------------------------------------
    # Write audit reports.
    # -------------------------------------------------------------

    with unmatched_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer_unmatched = csv.DictWriter(
            handle,
            fieldnames=[
                "csv_line",
                "label_hierarchy",
            ],
        )

        writer_unmatched.writeheader()
        writer_unmatched.writerows(
            unmatched_rows
        )

    with ambiguous_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer_ambiguous = csv.DictWriter(
            handle,
            fieldnames=[
                "csv_line",
                "label_hierarchy",
                "candidate_lookup_lines",
                "candidate_layers",
            ],
        )

        writer_ambiguous.writeheader()
        writer_ambiguous.writerows(
            ambiguous_rows
        )

    # -------------------------------------------------------------
    # Backup original and replace with enriched file.
    # -------------------------------------------------------------

    shutil.copy2(
        main_path,
        backup_path,
    )

    os.replace(
        temp_path,
        main_path,
    )

except Exception as error:
    if temp_path.exists():
        temp_path.unlink()

    raise SystemExit(
        f"ERROR: {error}"
    )


# =============================================================================
# 4. SUMMARY
# =============================================================================

print("")
print("VME TAXONOMY ENRICHMENT COMPLETE")
print("================================")
print(
    f"VME lookup rows loaded:     {len(lookup_rows):,}"
)
print(
    f"VME rows in BIIGLE file:    {vme_rows:,}"
)
print(
    f"Exact hierarchy matches:    {exact_matches:,}"
)
print(
    f"Safe fallback matches:      {fallback_matches:,}"
)
print(
    f"Unmatched VME rows:         {len(unmatched_rows):,}"
)
print(
    f"Ambiguous VME rows:         {len(ambiguous_rows):,}"
)
print(
    f"Non-VME rows untouched:     {non_vme_rows:,}"
)
print("")
print("Filled for matched VME rows:")
for column in TAXON_COLUMNS:
    print(f"  {column}")
print("")
print("Updated:")
print(f"  {main_path}")
print("")
print("Backup:")
print(f"  {backup_path}")
print("")
print("Audit reports:")
print(f"  {unmatched_path}")
print(f"  {ambiguous_path}")
print("")

PYTHON
