#!/usr/bin/env bash

set -euo pipefail

# =============================================================================
# FRAME-LEVEL SUBSTRATE COMPOSITION CLASS
# =============================================================================
#
# Run FROM the project folder containing:
#
#   combined_biigle_annotations.csv
#
# Recommended location:
#
#   ./scripts/add_frame_substrate_composition_class.sh
#
# Run:
#
#   chmod +x ./scripts/add_frame_substrate_composition_class.sh
#   ./scripts/add_frame_substrate_composition_class.sh
#
# PURPOSE
# -------
#
# Convert Abiotic POINT annotations into one useful whole-frame physical
# substrate predictor, while keeping relief as a separate frame predictor.
#
# Example:
#
#   soft substrate
#   62% silt
#   38% bioturbation
#
# becomes:
#
#   frame_substrate_class = soft_silt_bioturbation
#
# rather than:
#
#   mixed
#
# IMPORTANT
# ---------
#
# This script uses ONLY Abiotic POINT annotations.
#
# Biological annotations such as Macroalgae are NOT included in the substrate
# predictor. This avoids using part of the biological response as a predictor
# in later diversity/community models.
#
# Run this BEFORE splitting into:
#
#   point_annotations.csv
#   vme_annotations.csv
#
# because every row belonging to the same filename will inherit the same
# frame-level substrate and relief values.
#
# INPUT COLUMNS
# -------------
#
#   filename
#   annotation_id
#   top_level
#   shape_name
#   releif
#   substrate
#   type
#   size
#
# "releif" intentionally matches the current dataset spelling.
#
# SUBSTRATE CLASS LOGIC
# ---------------------
#
# 1. Work at unique filename + annotation_id level.
#
# 2. Determine the dominant broad substrate matrix:
#
#      soft
#      hard
#
# 3. Determine physical modifiers:
#
#      - if size is present, use size
#          silt, sand, pebble, cobble, boulder, bedrockreef, etc.
#
#      - otherwise use type, except "wentworth"
#          basalt, limestone, ice, shelldead, ripples,
#          bacterialmats, bioturbation, icescour, wall, etc.
#
# 4. Append up to TWO modifiers to the broad matrix.
#
#    A modifier must occur in at least 20% of the classifiable substrate
#    points in that frame to be included.
#
# Examples:
#
#      soft + 82% silt
#          -> soft_silt
#
#      soft + 62% silt + 38% bioturbation
#          -> soft_silt_bioturbation
#
#      soft + 55% silt + 30% boulder + 15% sand
#          -> soft_silt_boulder
#
#      hard + 48% boulder + 42% cobble + 10% pebble
#          -> hard_boulder_cobble
#
#      soft + no modifier reaching 20%
#          -> soft
#
# Relief is calculated independently. If two relief states are substantial,
# both can be retained:
#
#      flat
#      flat_moderate
#      moderate_high
#
# OUTPUT COLUMNS APPENDED TO MASTER CSV
# -------------------------------------
#
#   frame_substrate_class
#   frame_relief_class
#   frame_substrate_n
#   frame_substrate_dominance_pct
#
# A separate one-row-per-frame audit is also written:
#
#   frame_substrate_summary.csv
#
# BACKUP
# ------
#
#   combined_biigle_annotations.before_frame_substrate.csv
#
# =============================================================================


MAIN_CSV="combined_biigle_annotations.csv"

PROJECT_DIR="$(pwd)"
MAIN_PATH="${PROJECT_DIR}/${MAIN_CSV}"

if [[ ! -f "$MAIN_PATH" ]]; then
    echo "ERROR: Cannot find ${MAIN_CSV} in:" >&2
    echo "  ${PROJECT_DIR}" >&2
    echo "" >&2
    echo "Run this script from the project folder containing the CSV." >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required." >&2
    exit 1
fi


python3 - "$MAIN_PATH" <<'PYTHON'

import csv
import os
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path


main_path = Path(sys.argv[1]).resolve()

backup_path = main_path.with_name(
    "combined_biigle_annotations.before_frame_substrate.csv"
)

temp_path = main_path.with_name(
    "combined_biigle_annotations.frame_substrate_tmp.csv"
)

summary_path = main_path.with_name(
    "frame_substrate_summary.csv"
)


NEW_COLUMNS = [
    "frame_substrate_class",
    "frame_relief_class",
    "frame_substrate_n",
    "frame_substrate_dominance_pct",
]

MIN_COMPONENT_PROPORTION = 0.20
MAX_MODIFIERS = 2
MAX_RELIEF_COMPONENTS = 2


def clean_text(value):
    return str(value or "").replace("\ufeff", "").strip()


def norm(value):
    return clean_text(value).casefold()


def slug(value):
    value = norm(value)
    value = value.replace("&", "and")
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def find_index(headers, wanted, required=True):
    wanted_norm = norm(wanted)

    for index, header in enumerate(headers):
        if norm(header) == wanted_norm:
            return index

    if required:
        raise ValueError(
            f"Required column not found: {wanted}"
        )

    return None


def ranked_components(values):
    counts = Counter(
        value
        for value in values
        if value
    )

    return sorted(
        counts.items(),
        key=lambda item: (
            -item[1],
            item[0],
        ),
    )


def make_broad_matrix(values):
    """
    Return one broad substrate prefix.

    If there is a true tie for the most common broad substrate,
    retain both tied components rather than calling the frame "mixed".
    """
    ranked = ranked_components(values)

    if not ranked:
        return "", "", ""

    top_count = ranked[0][1]

    tied = [
        category
        for category, count in ranked
        if count == top_count
    ]

    prefix = "_".join(
        tied
    )

    dominance_pct = (
        100.0 * top_count / len(values)
    )

    return (
        prefix,
        " | ".join(tied),
        round(dominance_pct, 1),
    )


def make_modifier_list(values, total_substrate_n):
    """
    Keep up to two physical modifiers, in descending abundance order,
    where each accounts for >=20% of all classifiable substrate points.
    """
    if total_substrate_n <= 0:
        return [], []

    ranked = ranked_components(
        values
    )

    retained = []
    audit = []

    for category, count in ranked:
        proportion = (
            count / total_substrate_n
        )

        if proportion < MIN_COMPONENT_PROPORTION:
            continue

        retained.append(
            category
        )

        audit.append(
            (
                category,
                count,
                round(
                    100.0 * proportion,
                    1,
                ),
            )
        )

        if len(retained) >= MAX_MODIFIERS:
            break

    return retained, audit


def make_relief_class(values):
    """
    Keep up to two substantial relief states.

    Relief is independent of the substrate class.
    """
    values = [
        value
        for value in values
        if value
    ]

    if not values:
        return "", [], ""

    ranked = ranked_components(
        values
    )

    retained = []
    audit = []

    for category, count in ranked:
        proportion = (
            count / len(values)
        )

        if (
            proportion >= MIN_COMPONENT_PROPORTION
            or not retained
        ):
            retained.append(
                category
            )

            audit.append(
                (
                    category,
                    count,
                    round(
                        100.0 * proportion,
                        1,
                    ),
                )
            )

        if len(retained) >= MAX_RELIEF_COMPONENTS:
            break

    relief_class = "_".join(
        retained
    )

    primary_pct = (
        audit[0][2]
        if audit
        else ""
    )

    return (
        relief_class,
        audit,
        primary_pct,
    )


# =============================================================================
# 1. READ CSV
# =============================================================================

with main_path.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:

    reader = csv.reader(
        handle
    )

    try:
        header = next(
            reader
        )
    except StopIteration:
        raise SystemExit(
            "ERROR: combined_biigle_annotations.csv is empty."
        )

    rows = list(
        reader
    )


width = max(
    [len(header)]
    + [len(row) for row in rows]
)

if len(header) < width:
    header.extend(
        [""] * (
            width - len(header)
        )
    )

for row in rows:
    if len(row) < width:
        row.extend(
            [""] * (
                width - len(row)
            )
        )


# =============================================================================
# 2. REMOVE OLD GENERATED FRAME COLUMNS
# =============================================================================

generated_names = {
    norm(column)
    for column in NEW_COLUMNS
}

keep_indices = [
    index
    for index, name in enumerate(
        header
    )
    if norm(name) not in generated_names
]

header = [
    header[index]
    for index in keep_indices
]

rows = [
    [
        row[index]
        if index < len(row)
        else ""
        for index in keep_indices
    ]
    for row in rows
]


# Trim only trailing unnamed columns.
last_named_index = -1

for index, name in enumerate(
    header
):
    if clean_text(name):
        last_named_index = index

if last_named_index < 0:
    raise SystemExit(
        "ERROR: No named columns found."
    )

header = header[
    : last_named_index + 1
]

rows = [
    row[
        : last_named_index + 1
    ]
    for row in rows
]


# =============================================================================
# 3. FIND REQUIRED COLUMNS
# =============================================================================

filename_idx = find_index(
    header,
    "filename",
)

annotation_idx = find_index(
    header,
    "annotation_id",
)

top_level_idx = find_index(
    header,
    "top_level",
)

relief_idx = find_index(
    header,
    "releif",
)

substrate_idx = find_index(
    header,
    "substrate",
)

type_idx = find_index(
    header,
    "type",
)

size_idx = find_index(
    header,
    "size",
)

shape_idx = find_index(
    header,
    "shape_name",
    required=False,
)


# =============================================================================
# 4. RESOLVE UNIQUE ABIOTIC POINTS
# =============================================================================

point_records = defaultdict(
    lambda: {
        "broad_candidates": [],
        "modifier_candidates": [],
        "relief_candidates": [],
    }
)

frame_names = set()


for row_number, row in enumerate(
    rows,
    start=2,
):

    filename = clean_text(
        row[filename_idx]
    )

    if not filename:
        continue

    frame_names.add(
        filename
    )

    if norm(
        row[top_level_idx]
    ) != "abiotic":
        continue

    if (
        shape_idx is not None
        and clean_text(
            row[shape_idx]
        )
        and norm(
            row[shape_idx]
        ) != "point"
    ):
        continue

    annotation_id = clean_text(
        row[annotation_idx]
    )

    if not annotation_id:
        annotation_id = (
            f"__missing_annotation_id_row_{row_number}"
        )

    key = (
        filename,
        annotation_id,
    )

    broad = slug(
        row[substrate_idx]
    )

    substrate_type = slug(
        row[type_idx]
    )

    substrate_size = slug(
        row[size_idx]
    )

    relief = slug(
        row[relief_idx]
    )

    if broad:
        point_records[key][
            "broad_candidates"
        ].append(
            broad
        )

    if substrate_size:
        # Size is the finest physical modifier.
        point_records[key][
            "modifier_candidates"
        ].append(
            (
                3,
                substrate_size,
            )
        )

    elif (
        substrate_type
        and substrate_type != "wentworth"
    ):
        point_records[key][
            "modifier_candidates"
        ].append(
            (
                2,
                substrate_type,
            )
        )

    if relief:
        point_records[key][
            "relief_candidates"
        ].append(
            relief
        )


# =============================================================================
# 5. COLLAPSE MULTIPLE EXPORTED ROWS TO ONE OBSERVATION PER ANNOTATION_ID
# =============================================================================

frame_broad = defaultdict(
    list
)

frame_modifiers = defaultdict(
    list
)

frame_relief = defaultdict(
    list
)

ambiguous_broad_points = 0
ambiguous_modifier_points = 0
ambiguous_relief_points = 0


for (
    filename,
    annotation_id,
), record in point_records.items():

    broad_values = set(
        record[
            "broad_candidates"
        ]
    )

    if len(broad_values) == 1:
        frame_broad[
            filename
        ].append(
            next(
                iter(
                    broad_values
                )
            )
        )

    elif len(broad_values) > 1:
        ambiguous_broad_points += 1


    modifier_candidates = record[
        "modifier_candidates"
    ]

    if modifier_candidates:
        max_specificity = max(
            specificity
            for specificity, value
            in modifier_candidates
        )

        best_values = {
            value
            for specificity, value
            in modifier_candidates
            if specificity == max_specificity
        }

        if len(best_values) == 1:
            frame_modifiers[
                filename
            ].append(
                next(
                    iter(
                        best_values
                    )
                )
            )

        else:
            ambiguous_modifier_points += 1


    relief_values = set(
        record[
            "relief_candidates"
        ]
    )

    if len(relief_values) == 1:
        frame_relief[
            filename
        ].append(
            next(
                iter(
                    relief_values
                )
            )
        )

    elif len(relief_values) > 1:
        ambiguous_relief_points += 1


# =============================================================================
# 6. BUILD WHOLE-FRAME SUBSTRATE AND RELIEF CLASSES
# =============================================================================

frame_results = {}


for filename in sorted(
    frame_names
):

    broad_values = (
        frame_broad.get(
            filename,
            [],
        )
    )

    modifier_values = (
        frame_modifiers.get(
            filename,
            [],
        )
    )

    relief_values = (
        frame_relief.get(
            filename,
            [],
        )
    )

    substrate_n = len(
        broad_values
    )

    (
        broad_prefix,
        broad_primary,
        broad_dominance_pct,
    ) = make_broad_matrix(
        broad_values
    )

    (
        modifiers,
        modifier_audit,
    ) = make_modifier_list(
        modifier_values,
        substrate_n,
    )

    substrate_parts = []

    if broad_prefix:
        substrate_parts.append(
            broad_prefix
        )

    for modifier in modifiers:
        if (
            modifier
            and modifier not in substrate_parts
        ):
            substrate_parts.append(
                modifier
            )

    frame_substrate_class = (
        "_".join(
            substrate_parts
        )
        if substrate_parts
        else ""
    )

    (
        frame_relief_class,
        relief_audit,
        relief_primary_pct,
    ) = make_relief_class(
        relief_values
    )

    primary_modifier = (
        modifier_audit[0][0]
        if modifier_audit
        else ""
    )

    primary_modifier_pct = (
        modifier_audit[0][2]
        if modifier_audit
        else ""
    )

    secondary_modifier = (
        modifier_audit[1][0]
        if len(modifier_audit) > 1
        else ""
    )

    secondary_modifier_pct = (
        modifier_audit[1][2]
        if len(modifier_audit) > 1
        else ""
    )

    frame_results[
        filename
    ] = {
        "frame_substrate_class":
            frame_substrate_class,
        "frame_relief_class":
            frame_relief_class,
        "frame_substrate_n":
            str(
                substrate_n
            ),
        "frame_substrate_dominance_pct":
            (
                ""
                if primary_modifier_pct == ""
                else f"{primary_modifier_pct:.1f}"
            ),
        # Audit-only:
        "broad_primary":
            broad_primary,
        "broad_dominance_pct":
            broad_dominance_pct,
        "primary_modifier":
            primary_modifier,
        "primary_modifier_pct":
            primary_modifier_pct,
        "secondary_modifier":
            secondary_modifier,
        "secondary_modifier_pct":
            secondary_modifier_pct,
        "relief_primary_pct":
            relief_primary_pct,
    }


# =============================================================================
# 7. APPEND FRAME VARIABLES TO EVERY ROW
# =============================================================================

output_header = (
    header
    + NEW_COLUMNS
)

output_rows = []


for row in rows:

    filename = clean_text(
        row[filename_idx]
    )

    result = frame_results.get(
        filename,
        {},
    )

    output_rows.append(
        row
        + [
            result.get(
                "frame_substrate_class",
                "",
            ),
            result.get(
                "frame_relief_class",
                "",
            ),
            result.get(
                "frame_substrate_n",
                "0",
            ),
            result.get(
                "frame_substrate_dominance_pct",
                "",
            ),
        ]
    )


# =============================================================================
# 8. WRITE ONE-ROW-PER-FRAME AUDIT
# =============================================================================

summary_header = [
    "filename",
    "frame_substrate_class",
    "frame_relief_class",
    "frame_substrate_n",
    "frame_substrate_dominance_pct",
    "dominant_broad_substrate",
    "broad_substrate_dominance_pct",
    "primary_modifier",
    "primary_modifier_pct",
    "secondary_modifier",
    "secondary_modifier_pct",
    "relief_primary_pct",
]


with summary_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=summary_header,
    )

    writer.writeheader()

    for filename in sorted(
        frame_results
    ):

        result = frame_results[
            filename
        ]

        writer.writerow(
            {
                "filename":
                    filename,
                "frame_substrate_class":
                    result[
                        "frame_substrate_class"
                    ],
                "frame_relief_class":
                    result[
                        "frame_relief_class"
                    ],
                "frame_substrate_n":
                    result[
                        "frame_substrate_n"
                    ],
                "frame_substrate_dominance_pct":
                    result[
                        "frame_substrate_dominance_pct"
                    ],
                "dominant_broad_substrate":
                    result[
                        "broad_primary"
                    ],
                "broad_substrate_dominance_pct":
                    (
                        ""
                        if result[
                            "broad_dominance_pct"
                        ] == ""
                        else f"{result['broad_dominance_pct']:.1f}"
                    ),
                "primary_modifier":
                    result[
                        "primary_modifier"
                    ],
                "primary_modifier_pct":
                    (
                        ""
                        if result[
                            "primary_modifier_pct"
                        ] == ""
                        else f"{result['primary_modifier_pct']:.1f}"
                    ),
                "secondary_modifier":
                    result[
                        "secondary_modifier"
                    ],
                "secondary_modifier_pct":
                    (
                        ""
                        if result[
                            "secondary_modifier_pct"
                        ] == ""
                        else f"{result['secondary_modifier_pct']:.1f}"
                    ),
                "relief_primary_pct":
                    (
                        ""
                        if result[
                            "relief_primary_pct"
                        ] == ""
                        else f"{result['relief_primary_pct']:.1f}"
                    ),
            }
        )


# =============================================================================
# 9. WRITE AND VERIFY UPDATED MASTER CSV
# =============================================================================

with temp_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.writer(
        handle
    )

    writer.writerow(
        output_header
    )

    writer.writerows(
        output_rows
    )


with temp_path.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:

    checker = csv.reader(
        handle
    )

    check_header = next(
        checker
    )

    check_rows = sum(
        1
        for _ in checker
    )


if check_rows != len(rows):
    temp_path.unlink(
        missing_ok=True
    )

    raise SystemExit(
        "ERROR: Output row-count verification failed."
    )


if (
    check_header[
        -len(NEW_COLUMNS):
    ]
    != NEW_COLUMNS
):
    temp_path.unlink(
        missing_ok=True
    )

    raise SystemExit(
        "ERROR: New fields were not appended at the true end of the CSV."
    )


# =============================================================================
# 10. BACKUP + REPLACE
# =============================================================================

shutil.copy2(
    main_path,
    backup_path,
)

os.replace(
    temp_path,
    main_path,
)


# =============================================================================
# 11. REPORT
# =============================================================================

frames_with_substrate = sum(
    1
    for result in frame_results.values()
    if result[
        "frame_substrate_class"
    ]
)

print("")
print("FRAME SUBSTRATE COMPOSITION COMPLETE")
print("====================================")
print(f"Frames found:                   {len(frame_results):,}")
print(f"Frames with substrate class:    {frames_with_substrate:,}")
print(f"Ambiguous broad points skipped: {ambiguous_broad_points:,}")
print(f"Ambiguous modifiers skipped:    {ambiguous_modifier_points:,}")
print(f"Ambiguous relief points skipped:{ambiguous_relief_points:,}")
print("")
print("Appended:")
print("  frame_substrate_class")
print("  frame_relief_class")
print("  frame_substrate_n")
print("  frame_substrate_dominance_pct")
print("")
print("Updated:")
print(f"  {main_path}")
print("")
print("Backup:")
print(f"  {backup_path}")
print("")
print("Frame audit:")
print(f"  {summary_path}")
print("")

PYTHON
