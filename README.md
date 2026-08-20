# BIIGLE Benthic Image Analysis Pipeline

> [!NOTE]
> **Status: validated worked-example analytical pipeline.**  
> Processing Steps 00–07 and ecological Analyses 20–28 have been executed and statistically audited on the included worked example. Compact QC/result tables and GitHub-viewable PNG figures are versioned; large reconstructable datasets, distance matrices, model objects and duplicate PDF figures are regenerated locally and ignored by Git.

---

## Purpose

The aim of this repository is to make a complete BIIGLE benthic-image workflow **reproducible, inspectable and reusable by another researcher**. A user should be able to provide a manually concatenated BIIGLE annotation table plus a small set of lookup/metadata sheets, run numbered processing stages, inspect the output from every transformation, and then reproduce the ecological analyses and figures shown in this GitHub README.

The README is intended to function simultaneously as:

- a **data dictionary**;
- a **processing manual**;
- a **worked example** using the included dataset;
- a **statistical analysis guide**;
- a **figure gallery** linked directly to generated outputs; and
- a **template for adding future analyses** without changing upstream processing.

> [!IMPORTANT]
> The pipeline begins **after individual BIIGLE CSV exports have been manually concatenated**. Concatenation of the raw BIIGLE ZIP files is deliberately left as a manual prerequisite in the current workflow.

---

# Contents

1. [Getting Started on Mac](#getting-started-on-mac)
2. [Repository Structure](#repository-structure)
3. [Core Analytical Principles](#core-analytical-principles)
4. [Worked-Example Input Snapshot](#worked-example-input-snapshot)
5. [Required Input Files](#required-input-files)
6. [Processing Workflow](#processing-workflow)
7. [Final Datasets](#final-datasets)
8. [Running the Ecological Analyses](#running-the-ecological-analyses)
9. [Validated Worked-Example Results](#validated-worked-example-results)
10. [Statistical Analysis Map](#statistical-analysis-map)
11. [Alpha Diversity](#alpha-diversity)
12. [Community Resemblance and NMDS](#community-resemblance-and-nmds)
13. [PERMANOVA](#permanova)
14. [CAP / dbRDA](#cap--dbrda)
15. [Regression and Environmental Models](#regression-and-environmental-models)
16. [Hierarchical Univariate Group Comparisons](#hierarchical-univariate-group-comparisons)
17. [Taxon Responses](#taxon-responses)
18. [VME Analyses](#vme-analyses)
19. [Figure and GitHub-Linking System](#figure-and-github-linking-system)
20. [Central Configuration](#central-configuration)
21. [Quality Control](#quality-control)
22. [Adding a New Analysis](#adding-a-new-analysis)
23. [Reproducibility and Interpretation](#reproducibility-and-interpretation)
24. [Development Status](#development-status)

---

# Getting Started on Mac

The recommended workflow is to **work locally on a Mac inside a cloned GitHub repository** and use GitHub as the remote, versioned and shareable record of the project. Do not repeatedly download new copies of the repository into `Downloads/`. Maintain one local working copy and synchronise it with GitHub using Git.

Recommended local location:

```text
~/Documents/GitHub/BIIGLE-Benthic-Pipeline/
```

The complete installation, cloning, authentication, file-placement and day-to-day Git workflow is documented in:

[**SETUP.md — Getting Started on Mac**](SETUP.md)

A typical working session is:

```bash
cd ~/Documents/GitHub/BIIGLE-Benthic-Pipeline
git pull

# edit scripts, update Sheets/ when appropriate, and run analyses locally

git status

# Stage only the intended files; avoid blindly staging generated data.
git add Scripts/ README.md
git status
git commit -m "Describe the completed change"
git push
```

> [!IMPORTANT]
> `Raw/` is intended as **local provenance/archive storage** for original BIIGLE ZIP exports and other large source material. Its contents are ignored by Git by default. The reproducible pipeline begins from the manually concatenated `Sheets/combined_biigle_annotations.csv` plus the prerequisite metadata/lookup sheets.

> [!NOTE]
> The repository versions the canonical prerequisite `Sheets/`, compact worked-example QC/result tables, and PNG figures. `Data/Intermediate/`, generated `Data/Final/` tables, large distance matrices, model objects, row-level reconstructed analysis tables and duplicate PDF figures are reproducible products and are ignored by Git.

---

# Repository Structure

```text
BIIGLE-Benthic-Pipeline/
│
├── README.md
├── SETUP.md
├── .gitignore
│
├── Raw/
│   ├── README.md
│   └── original_BIIGLE_exports.zip   # local only; not committed
│
├── Sheets/
│   ├── combined_biigle_annotations.csv
│   ├── gl_latlong.csv
│   ├── depth_temp_frameid.csv
│   ├── taxon_list.csv
│   └── vme_taxon_list.csv
│
├── Scripts/
│   ├── 00_config.R
│   ├── 00_validate_inputs.sh
│   ├── 01_add_lat_long.sh
│   ├── 02_add_top_level.sh
│   ├── 03_add_taxonomy.sh
│   ├── 04_add_frame_metadata.sh
│   ├── 05_add_vme_taxonomy.sh
│   ├── 06_add_frame_substrate_class.sh
│   ├── 07_split_point_and_vme_data.sh
│   ├── 20_analysis_qc.R
│   ├── 21_analysis_alpha_diversity.R
│   ├── 22_analysis_nmds.R
│   ├── 23_analysis_permanova.R
│   ├── 24_analysis_cap_dbrda.R
│   ├── 25_analysis_environment.R
│   ├── 26_analysis_univariate.R
│   ├── 27_analysis_taxon_responses.R
│   └── 28_analysis_vme.R
│
├── Data/
│   ├── Intermediate/
│   │   ├── 01_annotations_latlong.csv
│   │   ├── 02_annotations_top_level.csv
│   │   ├── 03_annotations_taxonomy.csv
│   │   ├── 04_annotations_environment.csv
│   │   ├── 05_annotations_vme_taxonomy.csv
│   │   └── 06_annotations_frame_substrate.csv
│   │
│   └── Final/
│       ├── point_annotations.csv
│       └── vme_annotations.csv
│
├── Analyses/
│   ├── 00_Quality_Control/
│   ├── 01_Frame_Composition/
│   ├── 02_Alpha_Diversity/
│   │   ├── Shannon_Wiener/
│   │   ├── Simpson/
│   │   ├── Richness/
│   │   └── Rarefaction/
│   ├── 03_Community_Composition/
│   │   ├── Bray_Curtis/
│   │   └── NMDS/
│   ├── 04_PERMANOVA/
│   ├── 05_CAP_dbRDA/
│   ├── 06_Environmental_Associations/
│   ├── 07_Univariate_Tests/
│   ├── 08_Taxon_Responses/
│   └── 09_VME/
│
└── Figures/
    ├── Fig_01_...
    ├── Fig_02_...
    └── ...
```

| Folder | Meaning |
|---|---|
| `Raw/` | Local-only archive of original BIIGLE downloads; ignored by Git by default |
| `Sheets/` | Researcher-supplied prerequisite tables; treated as immutable |
| `Scripts/` | Numbered data-processing and analysis code |
| `Data/Intermediate/` | Stage-by-stage transformed datasets |
| `Data/Final/` | Stable analytical point and VME tables |
| `Analyses/` | Statistical tables, model objects, diagnostics and summaries |
| `Figures/` | Numbered plots linked directly into this README |

> [!IMPORTANT]
> **Never silently overwrite `Raw/` or `Sheets/`.** `Raw/` is a local provenance archive and `Sheets/` contains the canonical prerequisite inputs. A formal pipeline should read these inputs and write transformed files to `Data/Intermediate/`.

---

# Core Analytical Principles

## 1. The frame is the principal ecological sampling unit

The BIIGLE point design contains approximately **100 random annotation points per image frame**. Individual annotation points provide information used to characterise each frame, but inferential comparisons of biodiversity or assemblage structure should generally operate at the **frame level**, with the dependence among frames from the same dive/site explicitly considered.

## 2. Point and VME annotations are different sampling designs

`point_annotations.csv` represents the random-point sampling design used for cover, diversity and community analyses.

`vme_annotations.csv` contains VME annotations made when VME features were observed. These targeted records are **not equivalent random abundance samples** and therefore must not inflate point totals, percentage cover, richness or Bray–Curtis matrices.

## 3. A BIIGLE row is not automatically a unique point

A single `annotation_id` can appear in more than one exported row if multiple labels are attached to the same annotation. Therefore:

```text
same filename + same annotation_id = same annotation point
```

whereas:

```text
same filename + different annotation_id = different point
```

This rule is fundamental for point counts, diversity and community matrices.

## 4. Standardised biological unit

The preferred biological unit is resolved hierarchically as:

```text
common_id_full
    ↓ if missing
common_id_mid
    ↓ if missing
common_id_short
    ↓ if missing
label_name
```

These IDs are standardised annotation/morphotaxon units and should not automatically be described as formally identified species.

## 5. Physical substrate is a frame-level environmental predictor

Abiotic point annotations are summarised into a whole-frame `frame_substrate_class`, while relief is retained independently as `frame_relief_class`. The resulting variables describe the observed physical environment of the frame; they do **not** assert knowledge of the hidden substrate directly beneath every biological organism.

Biological habitat formers such as macroalgae should not be folded into the physical substrate predictor when the same biological point sample is used as the response.

---

# Worked-Example Input Snapshot

The worked example in this repository is built from a manually concatenated BIIGLE annotation table plus four researcher-supplied metadata/lookup sheets. The original BIIGLE annotation exports are retained locally in `Raw/` for provenance but are **not** required by the automated pipeline after concatenation.

| Input | Records | Primary role |
|---|---:|---|
| Raw BIIGLE annotation exports | 25 ZIP archives | Provenance for the manually concatenated example annotation table |
| `combined_biigle_annotations.csv` | 55,096 annotation rows | Manually concatenated BIIGLE annotations |
| BIIGLE frames | 525 unique filenames | Image-level sampling units |
| dives | 19 unique dive IDs | Spatial/sampling grouping |
| `gl_latlong.csv` | 28 lookup records | Dive/site coordinates |
| `depth_temp_frameid.csv` | 395 frame records | Frame depth, temperature and time |
| `taxon_list.csv` | 288 lookup records | Taxonomy, common IDs and abiotic descriptors |
| `vme_taxon_list.csv` | 58 lookup records | VME taxonomy and common IDs |

## Raw BIIGLE exports used for the worked example

The following 25 BIIGLE CSV image-annotation report archives form the raw provenance for the worked-example `combined_biigle_annotations.csv`:

```text
32643_csv_image_annotation_report.zip
32572_csv_image_annotation_report.zip
32574_csv_image_annotation_report.zip
32575_csv_image_annotation_report.zip
32581_csv_image_annotation_report.zip
32583_csv_image_annotation_report.zip
32584_csv_image_annotation_report.zip
32585_csv_image_annotation_report.zip
32586_csv_image_annotation_report.zip
32622_csv_image_annotation_report.zip
32625_csv_image_annotation_report.zip
32626_csv_image_annotation_report.zip
32627_csv_image_annotation_report.zip
32629_csv_image_annotation_report.zip
32630_csv_image_annotation_report.zip
32633_csv_image_annotation_report.zip
32634_csv_image_annotation_report.zip
32635_csv_image_annotation_report.zip
32636_csv_image_annotation_report.zip
32637_csv_image_annotation_report.zip
32638_csv_image_annotation_report.zip
32639_csv_image_annotation_report.zip
32640_csv_image_annotation_report.zip
32641_csv_image_annotation_report.zip
32642_csv_image_annotation_report.zip
```

These archives may be stored locally under `Raw/`. The current pipeline deliberately begins **after** these exports have been manually concatenated into `Sheets/combined_biigle_annotations.csv`; the ZIP files themselves are therefore not required to rerun the downstream enrichment and ecological analyses.

> [!IMPORTANT]
> The filenames shown throughout this repository are the **canonical local/GitHub filenames**. Browser or chat download suffixes such as `taxon_list(2).csv`, `depth_temp_frameid(4).csv`, or `combined_biigle_annotations(3).csv` are not part of the pipeline naming convention. Inside the repository these files must be named exactly `taxon_list.csv`, `depth_temp_frameid.csv`, `combined_biigle_annotations.csv`, and so on.

The raw combined table contains Biotic, Abiotic, Exclude, UNSURE, VME and occasional other annotation hierarchies. The formal processing stage derives `top_level` from the first component of `label_hierarchy`.

### Current metadata coverage in the worked example

- all 19 BIIGLE dive IDs match `gl_latlong.csv` after the documented case/hyphen/underscore normalisation;
- 209 of 525 unique BIIGLE frame filenames currently match `FrameID` values in `depth_temp_frameid.csv` exactly;
- unmatched frames are retained and reported explicitly rather than silently imputed.

> [!WARNING]
> Environmental metadata coverage is therefore incomplete in the current worked example. Analyses requiring depth or temperature must use only frames with valid matched metadata and must report the resulting sample size.

---

# Required Input Files

## `Sheets/combined_biigle_annotations.csv`

The manually concatenated BIIGLE export. The pipeline assumes the row structure produced by BIIGLE and retains `parent_csv` for provenance.

### Required raw columns

| Column | Type | Meaning |
|---|---|---|
| `dive_id` | character | Dive/sample identifier |
| `annotation_label_id` | ID | Exported annotation-label record |
| `label_id` | ID | BIIGLE label identifier |
| `label_name` | character | Human-readable annotation label |
| `label_hierarchy` | character | Full hierarchical annotation path |
| `user_id` | ID | BIIGLE user identifier |
| `firstname` | character | Annotator first name |
| `lastname` | character | Annotator surname |
| `image_id` | ID | BIIGLE image identifier |
| `filename` | character | Frame filename; principal frame join key |
| `image_longitude` | numeric/blank | BIIGLE image-level longitude if supplied |
| `image_latitude` | numeric/blank | BIIGLE image-level latitude if supplied |
| `shape_id` | ID | BIIGLE annotation shape type ID |
| `shape_name` | character | `Point`, `Polygon`, etc. |
| `points` | character | Pixel coordinates defining annotation geometry |
| `attributes` | JSON-like text | Image/annotation attributes |
| `annotation_id` | ID | **Unique annotation identity used to resolve duplicate exports** |
| `created_at` | datetime | Annotation creation timestamp |
| `parent_csv` | character | Original source CSV; provenance field |

---

## `Sheets/gl_latlong.csv`

Dive/site coordinate table.

| Column | Required | Meaning |
|---|:---:|---|
| `TripID` | No | Expedition/trip identifier |
| `Company` | No | Platform/operator |
| `Vessel` | No | Vessel name |
| `Country` | No | Country/region |
| `SampleID` | No | Numeric/local sample ID |
| `FullID` | **Yes** | Join key matched to `dive_id` after normalisation |
| `Date` | No | Sampling date |
| `Location` | No | Site name |
| `Lat` | **Yes** | Latitude |
| `Long` | **Yes** | Longitude |

### Join

```text
combined_biigle_annotations.csv : dive_id
                       ↕
gl_latlong.csv                  : FullID
```

The join should be case-insensitive and tolerant of underscore/hyphen differences. Conflicting duplicate coordinates must trigger an error rather than be guessed.

---

## `Sheets/depth_temp_frameid.csv`

Frame-level environmental/temporal metadata.

| Column | Required | Meaning |
|---|:---:|---|
| `id` | No | Source row/telemetry ID |
| `time` | No | Source timestamp |
| `Elapsed seconds` | No | Elapsed source-video time |
| `FrameID` | **Yes** | Exact frame filename join key |
| `video_time` | **Yes** | Video timestamp corresponding to frame |
| `depth_m` | **Yes** | Recorded depth; negative raw values are retained but analysed as absolute depth magnitude |
| `temperature` | **Yes** | Water temperature |
| `Pressabs(hPa)` | No | Pressure field |
| `lat` | No | Telemetry latitude when available |
| `lon` | No | Telemetry longitude when available |

### Join

```text
combined_biigle_annotations.csv : filename
                       ↕
depth_temp_frameid.csv          : FrameID
```

Unmatched frames remain in the dataset with blank environmental values and are reported in QC output.

---

## `Sheets/taxon_list.csv`

The principal biological and abiotic lookup table.

### Hierarchy fields

```text
Layer_01 ... Layer_08
```

### Biological fields

```text
cpc_codes
kingdom
phylum
class
order
family
taxonomic_resolution
morphology
common_id_short
common_id_mid
common_id_full
```

### Physical habitat fields

```text
releif
substrate
type
size
```

> [!NOTE]
> `releif` is intentionally retained with the spelling currently used by the source lookup table. Renaming it mid-pipeline would break reproducibility unless all scripts and data dictionaries are migrated together.

---

## `Sheets/vme_taxon_list.csv`

Lookup table for `top_level == VME` annotations.

### VME hierarchy

```text
Layer_01 ... Layer_05
```

### Enrichment fields

```text
cpc_codes
kingdom
phylum
class
order
family
taxonomic_resolution
morphology
common_id_short
common_id_mid
common_id_full
```

Matching should use the ordered VME hierarchy, with exact matching preferred and only conservative unambiguous fallback matching permitted.

---

# Processing Workflow

Each formal stage produces a **new numbered intermediate dataset** and leaves the researcher-supplied files in `Sheets/` unchanged. The canonical wrapper scripts now enforce this non-destructive contract. During the current development version, those wrappers call the already-tested implementations retained under `Scripts/legacy_current/`; this preserves validated matching logic while keeping the public workflow stable.

The entire processing phase can be run from the repository root with:

```bash
./Scripts/run_processing_pipeline.sh .
```

This executes Steps 00–07 in order and writes a processing manifest to:

```text
Analyses/00_Quality_Control/07_processing_pipeline_manifest.csv
```

## Step 00 — Validate inputs

**Script:** `Scripts/00_validate_inputs.sh`  
**Inputs:** all files in `Sheets/`  
**Outputs:** `Analyses/00_Quality_Control/`

Checks:

- required files exist;
- required columns exist;
- conflicting duplicate lookup keys;
- `dive_id ↔ FullID` coverage;
- `filename ↔ FrameID` coverage;
- frame-level environmental metadata coverage by dive;
- whether metadata coverage is complete, absent or partial within each dive;
- the total raw-export frame universe;
- the random-point-sampled frame universe;
- unique random-point `annotation_id` counts per frame;
- frames that do not meet the nominal 100-point design;
- VME-only frames with no random-point annotations;
- annotation IDs carrying conflicting top-level classifications;
- unmatched IDs are explicitly written to QC files.

For the current worked example, validation identifies two distinct analytical universes:

```text
Raw-export frame universe:         525 frames
Random-point-sampled universe:     522 frames
VME-only/no-random-point frames:     3 frames

Frames with depth/temp metadata:   209 / 525
Dives with complete metadata:        8 / 19
Dives with no frame metadata:       11 / 19
Partially matched dives:             0 / 19
```

Among the 522 point-sampled frames, 430 resolve to exactly 100 unique non-VME point annotations and 92 are flagged for review. These deviations are **reported rather than silently corrected**. Later community models use configurable minimum/maximum total-point thresholds and minimum assigned-biotic-point thresholds so severely incomplete or anomalously oversampled frames do not automatically enter inferential community analyses.

Run:

```bash
./Scripts/00_validate_inputs.sh .
```

## Step 01 — Add latitude/longitude

**Join:** `dive_id ↔ FullID`  
**Adds:** `lat`, `long`  
**Output:** `Data/Intermediate/01_annotations_latlong.csv`

## Step 02 — Derive `top_level`

For each row:

```text
Biotic > Echinoderm > ...  →  Biotic
Abiotic > Habitat Relief... →  Abiotic
VME > Coral > ...           →  VME
Exclude                      →  Exclude
```

**Output:** `Data/Intermediate/02_annotations_top_level.csv`

## Step 03 — Taxonomic and abiotic enrichment

Uses `taxon_list.csv` to populate standardised taxonomy, common IDs and physical substrate fields.

Rules:

- unsupported ranks remain blank rather than being invented;
- `Exclude` and `UNSURE` do not receive biological taxonomy;
- common IDs are carried forward exactly as defined by the lookup table.

**Output:** `Data/Intermediate/03_annotations_taxonomy.csv`

## Step 04 — Frame environmental metadata

**Join:** `filename ↔ FrameID`  
**Adds:** `video_time`, `depth_m`, `temperature`

**Output:** `Data/Intermediate/04_annotations_environment.csv`

## Step 05 — VME taxonomic enrichment

Only `top_level == VME` rows are enriched from `vme_taxon_list.csv`.

**Output:** `Data/Intermediate/05_annotations_vme_taxonomy.csv`

## Step 06 — Frame substrate characterisation

Abiotic **Point** annotations are resolved at the unique `filename + annotation_id` level and summarised into compact frame predictors.

Appended fields:

```text
frame_substrate_class
frame_relief_class
frame_substrate_n
frame_substrate_dominance_pct
```

Example:

```text
soft + silt + substantial bioturbation
        ↓
frame_substrate_class = soft_silt_bioturbation

frame_relief_class = flat
```

Substrate and relief remain separate so their ecological effects can be evaluated independently.

**Output:** `Data/Intermediate/06_annotations_frame_substrate.csv`

## Step 07 — Split analytical sampling designs

Creates:

```text
Data/Final/point_annotations.csv
Data/Final/vme_annotations.csv
```

The original/intermediate master table is retained.

## Current worked-example processing test

The complete Steps 00–07 processing workflow has been executed successfully on the example dataset. The current outputs are:

| Stage | Output | Rows | Columns |
|---|---|---:|---:|
| 01 | `Data/Intermediate/01_annotations_latlong.csv` | 55,096 | 21 |
| 02 | `Data/Intermediate/02_annotations_top_level.csv` | 55,096 | 22 |
| 03 | `Data/Intermediate/03_annotations_taxonomy.csv` | 55,096 | 37 |
| 04 | `Data/Intermediate/04_annotations_environment.csv` | 55,096 | 40 |
| 05 | `Data/Intermediate/05_annotations_vme_taxonomy.csv` | 55,096 | 40 |
| 06 | `Data/Intermediate/06_annotations_frame_substrate.csv` | 55,096 | 44 |
| 07a | `Data/Final/point_annotations.csv` | 51,955 | 44 |
| 07b | `Data/Final/vme_annotations.csv` | 3,141 | 44 |

Additional worked-example checks include:

- all 55,096 annotation rows receive latitude/longitude values;
- 209 unique frame filenames receive depth/temperature metadata;
- 3,141 VME rows are separated from the random-point table;
- VME hierarchy matching resolves 2,517 rows exactly and 622 using conservative unambiguous fallback matching, with two unmatched VME rows reported for review;
- 512 of 525 raw-export frames receive a frame substrate class;
- source files in `Sheets/` remain unchanged.

---

# Final Datasets

## `point_annotations.csv`

Used for random-point ecological analyses, including:

- point-class proportions;
- biotic cover;
- richness/diversity;
- frame community matrices;
- environmental modelling;
- Bray–Curtis resemblance;
- NMDS;
- PERMANOVA;
- CAP/dbRDA;
- taxon response analysis.

VME rows must not contribute additional random points.

## `vme_annotations.csv`

Used for:

- VME occurrence by frame/dive;
- VME-type summaries;
- VME relationships with depth, temperature, latitude/longitude and physical habitat;
- targeted VME visualisations.

Absence from `vme_annotations.csv` is interpreted as a frame-level non-detection **only when that frame belongs to the defined point/frame sampling universe and the VME annotation workflow was consistently applied**.

---

# Running the Ecological Analyses

After Steps 00–07 have produced the final point and VME datasets, install the required R packages once:

```bash
Rscript Scripts/install_R_packages.R
```

Then run the complete analysis suite:

```bash
./Scripts/run_analysis_pipeline.sh .
```

Or run individual modules:

```bash
Rscript Scripts/20_analysis_qc.R
Rscript Scripts/21_analysis_alpha_diversity.R
Rscript Scripts/22_analysis_nmds.R
Rscript Scripts/23_analysis_permanova.R
Rscript Scripts/24_analysis_cap_dbrda.R
Rscript Scripts/25_analysis_environment.R
Rscript Scripts/26_analysis_univariate.R
Rscript Scripts/27_analysis_taxon_responses.R
Rscript Scripts/28_analysis_vme.R
```

For a completely fresh rerun of processing **and** analyses:

```bash
./Scripts/run_full_pipeline.sh .
```

All user-editable analytical choices are centralised in `Scripts/00_config.R`. In particular, the current community workflow uses square-root-transformed relative abundance (`COMMUNITY_TRANSFORM = "sqrt_relative"`) so Bray–Curtis comparisons remain interpretable when some frames depart from the nominal 100-point effort. Community-model eligibility is also configurable using total-point and assigned-biotic-point thresholds.

Every statistical module writes:

- machine-readable result tables to its matching `Analyses/` subfolder;
- numbered figures to `Figures/`;
- sample sizes, model eligibility and/or numerical diagnostics where applicable;
- reconstructable model objects locally when useful.

Only compact audit/result tables and PNG figures are intended for routine Git versioning; large reconstructable products remain local and are governed by `.gitignore`.

---

# Validated Worked-Example Results

The tables below summarise the **validated worked example currently committed to the repository**. They are not universal ecological expectations and should not be transferred to another dataset without rerunning the pipeline.

> [!IMPORTANT]
> Several analyses use different analytical universes because the required metadata and sampling designs differ. Always report the sample size and grouping level for the specific analysis rather than referring to a single project-wide `n`.

## Processing and QC

The source annotation table contains **55,096 BIIGLE export rows** across **525 raw-export frames** and **19 dives**. Processing separates the records into **51,955 point-table rows** and **3,141 VME rows** without loss of source records.

For random-point sampling:

- **522 frames** contain at least one non-VME point annotation;
- **430/522 frames** contain exactly 100 unique non-VME annotation points;
- **92/522 frames** depart from the nominal 100-point design and are retained rather than silently forced to 100;
- **430 frames** meet the configured community/alpha-model eligibility thresholds;
- those 430 eligible community frames represent **18 dives**.

After resolving repeated BIIGLE export rows by `filename + annotation_id`, the unique non-VME point classes are:

| Point class | Unique annotation points |
|---|---:|
| Abiotic | 25,417 |
| Biotic | 14,810 |
| Exclude | 11,561 |
| Unsure | 139 |
| Conflict | 6 |

Of 14,810 unique biotic points, **14,808 resolve to one community unit**; two contain conflicting community-unit assignments and remain explicitly flagged.

Environmental coverage is structured rather than random: depth and temperature are available for **209 of 522 point-sampled frames**. Analyses using these predictors therefore report smaller sample sizes and never impute missing environmental values.

VME enrichment resolves 3,141 VME rows, with 2,517 exact hierarchy matches, 622 conservative unambiguous fallback matches and two unresolved records retained for review.

Relevant audit outputs:

- [`Analyses/00_Quality_Control/20_frame_qc_summary.csv`](Analyses/00_Quality_Control/20_frame_qc_summary.csv)
- [`Analyses/00_Quality_Control/20_point_class_counts.csv`](Analyses/00_Quality_Control/20_point_class_counts.csv)
- [`Analyses/00_Quality_Control/20_community_unit_resolution_summary.csv`](Analyses/00_Quality_Control/20_community_unit_resolution_summary.csv)
- [`Analyses/00_Quality_Control/20_frame_metadata_completeness.csv`](Analyses/00_Quality_Control/20_frame_metadata_completeness.csv)

![Point-count QC](Figures/Fig_01_Point_Count_QC.png)

![Frame composition by dive](Figures/Fig_02_Frame_Composition_By_Dive.png)

## Alpha diversity

All **522 point-sampled frames** remain in the master frame table. Of these:

- 513 contain resolved biota;
- 9 contain no resolved biota;
- 439 are eligible for rarefaction to 10 assigned biotic points;
- 430 meet the stricter alpha-model thresholds.

Median values among the relevant frames are:

| Metric | Median |
|---|---:|
| Assigned biotic points | 23 |
| Observed common-ID richness | 5 |
| Rarefied richness at 10 biotic points | 3.771 |
| Shannon diversity | 1.162 |
| Inverse Simpson diversity | 2.572 |
| Pielou evenness | 0.785 |

Zero-biota frames are retained with observed richness = 0; diversity indices that are undefined without biota remain `NA`.

![Observed richness versus depth](Figures/Fig_03_Richness_Depth.png)

![Shannon diversity versus depth](Figures/Fig_04_Shannon_Depth.png)

![Shannon diversity by substrate](Figures/Fig_05_Shannon_Substrate.png)

## Community composition and NMDS

The main community matrix contains **430 eligible frames** and **73 resolved biological groups**, using square-root-transformed relative point abundance before Bray–Curtis dissimilarity.

The two-dimensional NMDS solution has stress **0.212**, whereas the three-dimensional solution has stress **0.158**, a **25.8% reduction**. The 3-D solution is therefore the more defensible low-dimensional representation; 2-D figures are retained as simplified visualisations and should be reported with their stress value.

Exploratory `envfit` results identify:

- temperature: `R² = 0.357`, BH-adjusted `p = 0.0016` across the frame-level environmental-fit family;
- substrate class: `R² = 0.492`, BH-adjusted `p = 0.0236`;
- depth: BH-adjusted `p = 0.103`;
- relief: BH-adjusted `p = 0.680`;
- latitude and longitude: no supported dive-level association.

These `envfit` results describe association with the ordination configuration; they are not interpreted as causal environmental drivers.

![NMDS latitude visualisation](Figures/Fig_06_NMDS_Latitude.png)

![NMDS by substrate](Figures/Fig_07_NMDS_Substrate.png)

![NMDS by relief](Figures/Fig_08_NMDS_Relief.png)

## PERMANOVA and PERMDISP

Frame-level PERMANOVA permutations are restricted by `dive_id`; latitude and longitude are tested after aggregation to the dive level.

Standalone blocked PERMANOVA gives the strongest raw signals for:

| Predictor | Frames | Dives | R² | Raw p | Global BH q |
|---|---:|---:|---:|---:|---:|
| Temperature | 179 | 7 | 0.143 | 0.0138 | 0.0517 |
| Substrate | 366 | 17 | 0.136 | 0.0155 | 0.0517 |
| Depth | 179 | 7 | 0.070 | 0.4558 | 0.6614 |
| Relief | 420 | 17 | 0.025 | 0.7920 | 0.7920 |

Under the predeclared **global BH correction across reported inferential PERMANOVA terms**, no standalone term crosses `q < 0.05`.

Categorical PERMANOVA must be interpreted with PERMDISP:

- substrate dispersion differs strongly among groups (`q = 0.0001`);
- relief dispersion also differs (`q = 0.0001`).

Therefore the substrate PERMANOVA signal cannot be attributed solely to different group centroids; heterogeneous multivariate dispersion may contribute. Relief shows no centroid-composition support despite strong dispersion differences.

Dive-level latitude and longitude PERMANOVAs are unsupported.

## CAP / dbRDA

The complete-case constrained community model contains **141 frames from 6 dives** and uses depth, temperature, substrate and relief together.

The model explains:

- raw constrained `R² = 0.381`;
- adjusted `R² = 0.353`;
- overall dive-blocked permutation `p = 0.038`.

Collinearity diagnostics are acceptable (`max VIF = 6.80`, zero VIFs > 10, zero aliased constraints).

Marginal term tests within this four-predictor dbRDA are:

| Predictor | p | BH q within dbRDA term family |
|---|---:|---:|
| Temperature | 0.0093 | **0.0372** |
| Substrate | 0.0256 | 0.0512 |
| Relief | 0.1934 | 0.2579 |
| Depth | 0.6964 | 0.6964 |

Temperature is therefore the most consistently supported measured environmental correlate of community composition in this worked example. The complete-case model is nevertheless based on only six dives and should not be presented as explaining 35% of assemblage variation across the entire survey.

The dive-level latitude + longitude dbRDA has raw `R² = 0.108` but adjusted `R² = -0.011`; neither coordinate is supported.

![Environmental dbRDA](Figures/Fig_09_CAP_dbRDA_Environment.png)

![Dive-level geographic dbRDA](Figures/Fig_10_CAP_dbRDA_Dive_Latitude.png)

## Alpha-diversity responses to continuous environment

Depth and temperature are analysed at the frame level with **within-dive centred predictors**. Linear/log-linear/GAM models include a dive random intercept, and Spearman significance is obtained by within-dive permutation. Latitude and longitude are analysed with the dive as the inferential unit.

Temperature shows a coherent negative association with local alpha diversity. After BH correction within the blocked-Spearman family:

| Response | Spearman ρ | BH q |
|---|---:|---:|
| Observed richness | -0.192 | 0.0120 |
| Rarefied richness | -0.212 | 0.0030 |
| Shannon diversity | -0.219 | 0.0030 |
| Inverse Simpson diversity | -0.175 | 0.0165 |

Hierarchical GAMs support the same broad pattern for observed richness, rarefied richness, Shannon and inverse Simpson diversity (`q ≈ 0.013–0.022` within the GAM family). The smooth effective degrees of freedom are modest, suggesting mild non-linearity rather than highly complex response shapes.

Model `R²`/deviance-explained values include the **dive random effect as well as the environmental term** and must not be described as variance explained by temperature alone.

Depth provides much weaker evidence for ordinary alpha diversity. One exploratory GAM association is retained for overall biotic point percentage versus depth (`q = 0.0394`); its ecological direction should be interpreted from the fitted response rather than from model fit alone.

Latitude and longitude show no consistent dive-level diversity signal.

## Hierarchical substrate and relief comparisons

Ordinary Welch tests/ANOVA are **not used inferentially** for repeated frame observations because frames are nested within dives.

Instead, Analysis 26 uses partial-F tests after controlling for dive, with factor labels permuted **within dive**. Pairwise substrate tests are restricted further to dives where both compared substrate classes co-occur.

Substrate is eligible for hierarchical inference, but none of the tested alpha-diversity/biotic-cover responses is supported after BH correction. Partial `R²` values are small (approximately 0.6–2.1%).

Relief is **not inferentially eligible** for these alpha-diversity comparisons because only one eligible dive contains within-dive relief variation.

This is an important distinction from the multivariate result: substrate can be related to **which taxa comprise an assemblage and to assemblage heterogeneity** without changing local richness or Shannon diversity.

## Focal-taxon responses

Analysis 27 selects the 12 most prevalent focal community units from the **430 QC-eligible community frames**. The most prevalent include:

| Focal community unit | Frame prevalence |
|---|---:|
| Brittle star | 90.7% |
| Burrowing bivalve | 58.8% |
| Branching hydrocoral | 44.2% |
| Annelid tube | 36.5% |
| Encrusting demosponge | 27.9% |

Depth and temperature models use within-dive centred predictors and a dive random intercept. Occurrence and proportional valid-point cover are analysed separately. Model warnings and convergence are recorded per fit; the validated rerun produced **no numerically unstable models requiring exclusion**.

The strongest corrected proportional-cover responses are:

- brittle-star cover increases with depth (`q = 0.00128`);
- brittle-star cover decreases with temperature (`q = 0.00100`);
- burrowing-clam cover increases with temperature (`q = 0.0186`).

This indicates taxon turnover in relative dominance rather than a uniform response of all biological groups to temperature.

No focal-taxon substrate response survives its declared BH correction family. Relief is not eligible for focal-taxon factor inference under the available within-dive replication. No focal taxon shows a supported latitude/longitude prevalence association after correction.

![Focal-taxon prevalence](Figures/Fig_11_Focal_Taxon_Prevalence.png)

## VME occurrence

VME analyses remain separate from the random-point community design.

Across the **522 screened point-frame universe**, 371 frames contain at least one VME annotation (frame prevalence **71.1%**). VME frame prevalence varies strongly among dives.

The most widespread VME units include hydrocorals, brittle stars, encrusting demosponges and branching hydrocorals. Raw VME annotation counts are retained for provenance but are **not interpreted as organism abundance**, because VME annotations are targeted features rather than random points.

For frames with matched depth/temperature metadata, the VME occurrence models contain **209 frames from 8 dives**:

| Predictor | Estimate | p | BH q | Direction |
|---|---:|---:|---:|---|
| Depth | +0.1108 | 0.00157 | **0.00314** | VME occurrence increases with depth |
| Temperature | -1.4211 | 0.02754 | **0.02754** | VME occurrence decreases with temperature |

These are separate models, so they are interpreted as associations rather than independent effects conditional on one another.

Substrate is unsupported for VME occurrence. Relief is supported within the VME analytical universe (`q = 0.0139`), which contains sufficient within-dive relief variation for this response even though relief was ineligible in the stricter alpha-diversity universe.

Latitude and longitude are unsupported at the dive level. The validated VME run produced no numerically unstable models requiring exclusion.

> [!CAUTION]
> Treating a frame without a VME annotation as a true VME absence is valid only if VME screening/annotation was applied consistently across the analysed frame universe. If screening was opportunistic, VME occurrence models must be downgraded to descriptive summaries.

![VME occurrence by dive](Figures/Fig_12_VME_Occurrence_By_Dive.png)

## Overall worked-example interpretation

Across the validated modules, the most coherent result is that **fine-scale environmental conditions are more informative than broad latitude/longitude for the observed benthic assemblages**. Temperature is the most consistent measured environmental correlate of whole-community composition and local diversity, while focal taxa show different response directions, indicating reorganisation rather than a uniform assemblage-wide decline. Substrate is associated more clearly with multivariate community structure/heterogeneity than with local alpha diversity. VME occurrence is widespread and is associated with depth, temperature and relief within its own sampling universe.

These are observational associations, not causal effects. Metadata incompleteness, the limited number of independent dives in some complete-case models, and the distinction between random-point and targeted-VME sampling must remain explicit in any manuscript interpretation.

---

# Statistical Analysis Map

| Analysis | Primary response | Predictors / grouping | Main purpose |
|---|---|---|---|
| Frame QC | point counts, Exclude %, metadata completeness | frame, dive, depth/temp | Verify sampling and data quality |
| Observed richness | number of common IDs | environment | Quantify local biological richness |
| Rarefied richness | richness standardised to equal biotic-point effort | environment | Compare richness at standard effort |
| Shannon diversity | abundance + richness/evenness | environment | Compare diversity structure |
| Inverse Simpson | dominance-weighted diversity | environment | Detect changes driven by common groups |
| Pielou evenness | evenness | environment | Assess distribution of points among groups |
| Hierarchical environmental models | diversity/cover | depth/temp within dive; lat/long among dives | Estimate continuous environmental associations without frame pseudoreplication |
| Hierarchical factor tests | diversity/cover | substrate/relief within dive | Test categorical responses while controlling for dive |
| Bray–Curtis | multivariate abundance matrix | pairwise frame resemblance | Quantify assemblage dissimilarity |
| NMDS | Bray–Curtis matrix | visual overlays/environment | Visualise community turnover |
| PERMANOVA | Bray–Curtis matrix | depth/temp/substrate/relief; dive-level geography | Test multivariate assemblage differences |
| PERMDISP | distance to multivariate median | categorical predictor | Diagnose dispersion differences relevant to PERMANOVA |
| CAP/dbRDA | community resemblance | constrained environment | Quantify/test community variation explained by environment |
| Focal-taxon models | occurrence and valid-point cover | environment with dive hierarchy | Identify taxa associated with gradients |
| VME occurrence | VME frame presence | environment with dive hierarchy | Evaluate VME distribution separately from point abundance |

> [!IMPORTANT]
> The repository makes multiple analytical approaches available, but it does not blindly run every mathematically possible comparison. A test is inferential only when the response/predictor types, replication, sampling hierarchy and numerical diagnostics make the comparison defensible. Ineligible and unstable models are recorded explicitly rather than silently omitted.

# Alpha Diversity

## Observed common-ID richness

### What are we testing?

How many distinct standardised biological units occur among the valid biological point observations within each frame?

### Response

```text
common_id_richness
```

### Unit

One frame.

### Interpretation

Higher values indicate more observed biological groups, but raw richness is sensitive to the number of biological points available. Therefore observed richness should be interpreted alongside sampling effort and rarefied richness.

---

## Rarefied richness

### What are we testing?

Whether frames differ in richness after standardising comparison to the same number of assigned biological points.

### Response

```text
rarefied_common_id_richness
```

### Why?

A frame with 60 biological points has more opportunities to encounter groups than a frame with 10. Rarefaction reduces this unequal-effort effect.

### Interpretation

A higher rarefied value means greater expected richness at the chosen common sampling effort. It does not estimate total unseen species richness.

---

## Shannon–Wiener diversity

### Question

Do frames differ in diversity when both **number of groups and distribution of observations among groups** are considered?

For relative proportions \(p_i\):

$$
H'=-\sum_i p_i\ln(p_i)
$$

### Response

```text
shannon_common_id_diversity
```

### Interpretation

- larger `H′` = a richer and/or more even assemblage;
- lower `H′` = fewer groups and/or stronger dominance;
- a significant environmental relationship means Shannon diversity changes along that predictor, not necessarily that species richness alone changes.

---

## Inverse Simpson diversity

### Question

Does diversity change when relatively common/dominant biological units receive greater weight?

### Interpretation

Inverse Simpson is less sensitive than Shannon to rare groups. Higher values indicate greater effective diversity and lower dominance by one/few groups.

---

## Pielou evenness

### Question

How evenly are biological observations distributed among the groups present?

$$
J'=\frac{H'}{\ln(S)}
$$

where `S` is observed richness.

Values closer to 1 indicate greater evenness.

---

# Community Resemblance and NMDS

## Community matrix

The default frame × biological-unit matrix is built from **unique resolved biological point annotations**.

```text
rows    = frames
columns = common_id units
cells   = point counts
```

### Default transformation

The current worked-example default is:

```text
frame counts
    ↓ convert to within-frame relative abundance
relative abundance
    ↓ square-root transform
sqrt(relative abundance)
```

This is configured as:

```r
COMMUNITY_TRANSFORM <- "sqrt_relative"
```

The relative-abundance step reduces sensitivity to unequal total point effort among the minority of frames that depart from the nominal 100-point design, while the square-root step reduces the influence of highly dominant biological groups. Alternative transformations (`none`, `sqrt`, `relative`, `sqrt_relative`) remain configurable in `Scripts/00_config.R`.

## Bray–Curtis resemblance

### What are we testing/quantifying?

Bray–Curtis measures compositional dissimilarity between every pair of frames based on their biological abundances.

Values range approximately from:

```text
0 = identical composition
1 = no shared abundance structure
```

### Why Bray–Curtis?

It is widely used for ecological abundance data and does not treat joint absences as evidence of similarity.

### Output

```text
Analyses/03_Community_Composition/Bray_Curtis/
```

including the matrix itself and any diagnostics needed downstream.

---

## NMDS — Non-metric Multidimensional Scaling

### Ecological question

**How does whole-community composition vary among frames, dives and environmental gradients?**

### Input

1. frame × common-ID abundance matrix;
2. optional square-root transformation;
3. Bray–Curtis dissimilarity.

### What NMDS does

NMDS attempts to arrange frames in a low-dimensional space while preserving the **rank order** of ecological dissimilarities. Frames plotted close together have more similar assemblages; frames far apart have less similar assemblages.

### NMDS stress

Stress measures how faithfully the reduced-dimensional ordination represents the dissimilarity ranks.

As a practical guide only:

| Stress | Broad interpretation |
|---:|---|
| `< 0.05` | excellent representation |
| `0.05–0.10` | very good |
| `0.10–0.20` | potentially useful with caution |
| `> 0.20` | weak representation; consider more dimensions or cautious interpretation |

Stress is not a p-value.

### Default visualisation requested for this pipeline

Points can be coloured by latitude using a continuous gradient:

```text
lower/southern latitude → lighter
higher/northern latitude → darker
```

This is controlled centrally rather than hard-coded into the NMDS script.

Example configuration:

```r
COLOUR_MODE <- "palette"
LATITUDE_PALETTE <- "viridis"
LATITUDE_DIRECTION <- -1
```

or custom endpoints:

```r
COLOUR_MODE <- "custom"
LATITUDE_SOUTH_COLOUR <- "#..."
LATITUDE_NORTH_COLOUR <- "#..."
```

Other overlays may include:

- depth;
- temperature;
- dive/site;
- `frame_substrate_class`;
- `frame_relief_class`;
- latitude/longitude;
- VME presence where analytically appropriate.

### Environmental vectors

`vegan::envfit()` may be used to fit continuous environmental gradients or categorical centroids onto an unconstrained NMDS.

**Interpretation:** a longer vector indicates a stronger ordination-space correlation; direction indicates where values increase. Permutation p-values test whether the fitted association is stronger than expected under the permutation scheme.

> [!CAUTION]
> NMDS axes themselves do not have fixed ecological meanings. Do not interpret “NMDS1” as depth unless an environmental fit/model supports that relationship.

### Planned GitHub figure

```markdown
![NMDS coloured by latitude](Figures/Fig_XX_NMDS_Latitude.png)
```

---

# PERMANOVA

## What are we testing?

PERMANOVA tests whether multivariate community composition differs in relation to one or more predictors.

Default response:

```text
Bray–Curtis dissimilarity among frame-level biological assemblages
```

The pipeline deliberately separates **frame-level** and **dive-level** predictors.

### Frame-level predictors

```text
depth_m
temperature
frame_substrate_class
frame_relief_class
```

These are tested at the frame level, with permutations restricted within `dive_id` where the design supports it.

### Dive-level predictors

```text
lat
long
```

Latitude and longitude are generally constant for all frames belonging to the same dive in this example. They therefore cannot be meaningfully tested as frame-level predictors while simultaneously restricting permutations within dive. For latitude/longitude PERMANOVA, the community matrix is first aggregated to the **dive level**, giving one assemblage sample per dive.

### Example hypotheses

**Depth:** community composition does not vary systematically with depth within the sampled dive structure.

**Temperature:** community composition does not vary systematically with temperature within the sampled dive structure.

**Substrate:** centroids of community composition do not differ among frame substrate classes.

**Dive-level latitude:** aggregated dive assemblage composition does not vary systematically with latitude.

**Combined frame model:** after accounting for the other included frame-level terms, the tested predictor explains no additional assemblage variation.

### Example frame-level model

```r
adonis2(
  community ~ depth_m + temperature +
    frame_substrate_class + frame_relief_class,
  data = frame_metadata,
  permutations = 9999,
  by = "margin",
  strata = frame_metadata$dive_id
)
```

Latitude and longitude are handled separately at the dive level rather than being included in this within-dive permutation model. The exact formula and permutation count are configurable.

### Outputs to interpret

| Output | Meaning |
|---|---|
| pseudo-F | separation attributable to predictor relative to residual variation |
| `R²` | proportion of multivariate variation attributed to the term/model |
| p-value | permutation evidence against the null model |

### Interpretation example

If substrate returns:

```text
R² = 0.12
p = 0.002
```

this means substrate class is associated with approximately 12% of the variation represented by the tested multivariate model and the observed association was uncommon under the permutation null. It does **not** prove that substrate causally determines the assemblage.

## PERMDISP must accompany categorical PERMANOVA

PERMANOVA can be sensitive to differences in within-group dispersion. Therefore categorical tests should be accompanied by `betadisper()` / permutation testing.

### Question

Are groups different because their **centroids** differ, or because one group is simply more variable/dispersed?

A significant PERMDISP result means the PERMANOVA must be interpreted cautiously because heterogeneity of multivariate dispersion may contribute to the result.

## Repeated frames / permutation restriction

Frames from the same dive are not automatically independent. The default design therefore provides a permutation-strata option using `dive_id`. Where the biological question requires inference beyond sampled dives, mixed/hierarchical models or alternative sampling-unit definitions may be necessary.

---

# CAP / dbRDA

CAP here refers to constrained analysis of community resemblance. The current pipeline implements the worked example using `vegan::dbrda()` with Bray–Curtis dissimilarity and reports overall, term-level and constrained-axis permutation tests.

## What are we testing?

**How much of the multivariate community pattern can be explained by specified environmental variables, and along what constrained gradients do assemblages separate?**

### Response

Bray–Curtis community resemblance.

### Candidate explanatory variables

The worked-example frame-level constrained ordination uses:

```text
depth_m
temperature
frame_substrate_class
frame_relief_class
```

Latitude/longitude are not placed into the same frame-level constrained model because they are primarily dive-level variables in this example dataset.

### Example

```r
dbrda(
  community ~ depth_m + temperature + frame_substrate_class +
    frame_relief_class,
  data = frame_metadata,
  distance = "bray"
)
```

### Tests

Permutation tests should be produced for:

- overall constrained model;
- individual terms;
- constrained axes where useful.

### Interpretation

Unlike NMDS, CAP/dbRDA is **constrained by the supplied environmental predictors**. Separation along a constrained axis therefore represents community variation associated with the environmental model, not merely the strongest unconstrained structure in the community matrix.

Use CAP alongside, not instead of, NMDS:

```text
NMDS = What compositional structure exists?
CAP  = What part of that structure is associated with specified predictors?
```

---

# Regression and Environmental Models

## Spearman rank correlation

### Question

Is there a monotonic relationship between two continuous variables without assuming a linear normal relationship?

Examples:

```text
Shannon diversity vs depth
richness vs temperature
biotic cover vs depth
focal-group abundance vs latitude
```

### Outputs

- `rho`: strength/direction of monotonic association;
- p-value;
- sample size.

A positive `rho` indicates increasing ranks together; a negative value indicates opposite directions.

---

## Linear regression

### Question

How does a continuous ecological response change with one or more predictors under an approximately linear relationship?

Example:

```r
lm(shannon_common_id_diversity ~ depth_m + temperature,
   data = frame_data)
```

### Interpret

- coefficient sign = direction;
- coefficient magnitude = expected change in response per predictor unit, conditional on other terms;
- confidence interval = uncertainty;
- `R²` = proportion of response variance explained by the model;
- p-value = evidence against a zero coefficient under model assumptions.

Diagnostics should assess residual patterns, influential observations and variance behaviour.

---

## Log-transformed linear regression

A log transformation can be useful for positive, right-skewed ecological responses or multiplicative relationships.

Example:

```r
lm(log1p(response) ~ depth_m, data = frame_data)
```

### Interpretation

Coefficients occur on the transformed scale. The README/output should state the exact transformation (`log`, `log10`, `log1p`, etc.) and, where possible, provide back-transformed predictions for biological interpretation.

> [!WARNING]
> Transformations are chosen because they improve model suitability or answer a specific multiplicative question—not simply because an untransformed p-value was non-significant.

---

## GAM — Generalised Additive Model

### Question

Does an ecological response vary **nonlinearly** along depth, temperature or another continuous environmental gradient?

Example:

```r
gam(
  shannon_common_id_diversity ~ s(depth_m, k = 4),
  data = frame_data,
  method = "REML"
)
```

### Interpretation

The smooth shows the estimated shape of the relationship. Effective degrees of freedom (`edf`) indicate complexity: values near 1 approximate a linear effect, while larger values indicate increasing nonlinearity.

Plots should show the fitted trend and uncertainty rather than relying only on p-values.

---

# Hierarchical Univariate Group Comparisons

Analysis 26 tests whether frame-level alpha diversity or biological cover differs among categorical habitat groups **without treating frames from the same dive as independent replicates**.

For substrate and relief:

1. factor levels must meet the configured minimum frame count;
2. at least two dives must be represented;
3. at least two dives must contain within-dive variation in the factor;
4. the response is modelled after controlling for dive as a blocking factor;
5. significance is obtained by permuting factor labels within dive.

Pairwise substrate comparisons are stricter: the two factor levels must co-occur within at least two dives, and minimum replication must still be met after restricting to those co-occurring dives.

`dive_id` itself is summarised descriptively rather than tested with individual frames as if they were independent dive replicates.

> [!IMPORTANT]
> Ordinary Welch t-tests and classical frame-level ANOVA are therefore **not inferentially eligible for the current worked-example sampling design**. The hierarchical blocked procedure is the implemented default.

Multiple global and pairwise tests use Benjamini–Hochberg false-discovery-rate control according to the result families recorded in the output tables.

Relevant outputs:

- `Analyses/07_Univariate_Tests/26_factor_test_eligibility.csv`
- `Analyses/07_Univariate_Tests/26_global_hierarchical_factor_tests.csv`
- `Analyses/07_Univariate_Tests/26_pairwise_test_eligibility.csv`
- `Analyses/07_Univariate_Tests/26_pairwise_hierarchical_factor_tests.csv`
- `Analyses/07_Univariate_Tests/26_dive_descriptive_summaries.csv`

# Taxon Responses

Analysis 27 asks which focal community units change in **frame occurrence** or **proportional valid-point cover** along environmental gradients.

Focal taxa are selected from the configured number of most prevalent resolved community units among QC-eligible community frames.

## Continuous predictors

For depth and temperature:

- the sampling unit is the frame;
- predictors are centred within dive;
- models include `s(dive_id, bs = "re")` as a dive random intercept;
- binary frame occurrence is modelled with binomial GAMs;
- proportional valid-point cover is modelled with quasibinomial models;
- linear and smooth forms are retained only when the fit is numerically stable.

Model warnings and convergence state are captured per fit. A critical convergence warning triggers a retry with an alternative optimizer; a model that remains unstable is written to an exclusion table and is not treated as inferential.

## Categorical predictors

Substrate and relief models are restricted to dives where the factor varies within dive and to sufficiently replicated factor levels. Rank-deficient tests are flagged and not interpreted.

## Geography

Latitude and longitude are effectively dive-level predictors. Focal-taxon frame prevalence is therefore aggregated to one value per dive before geographic association tests.

Relevant audit outputs include:

- `27_continuous_model_eligibility.csv`
- `27_continuous_taxon_models.csv`
- `27_factor_model_eligibility.csv`
- `27_factor_taxon_models.csv`
- `27_geography_taxon_associations.csv`
- `27_model_fit_diagnostics.csv`
- `27_models_excluded_numerical_instability.csv`

# VME Analyses

VME annotations are analysed separately from the random-point abundance matrix because they are **targeted annotations**, not additional random sampling points.

The primary inferential response is frame-level VME occurrence:

```text
has_vme = 1 if one or more VME annotations occur in the screened frame
has_vme = 0 otherwise
```

Raw VME polygon/annotation counts are retained for provenance and descriptive summaries but are not automatically interpreted as abundance.

## Depth and temperature

Continuous VME occurrence models:

- use the frame as the sampling unit;
- centre the environmental predictor within dive;
- include a dive random intercept;
- capture convergence warnings explicitly.

Depth and temperature are currently fitted in separate models. Their coefficients therefore describe associations, not mutually adjusted independent effects.

## Substrate and relief

Categorical VME models are restricted to dives where the factor varies within dive. Rank-deficient models are flagged rather than interpreted.

## Latitude and longitude

VME frame prevalence is aggregated to one value per dive before geographic tests.

> [!CAUTION]
> A frame with no recorded VME can be treated as a VME absence only when VME screening/annotation was applied consistently across the analysed frame universe. If VME annotation was opportunistic, occurrence models should be treated as descriptive rather than inferential.

Relevant outputs include:

- `28_vme_presence_by_dive.csv`
- `28_vme_type_prevalence.csv`
- `28_continuous_model_eligibility.csv`
- `28_continuous_vme_models.csv`
- `28_factor_model_eligibility.csv`
- `28_factor_vme_models.csv`
- `28_geography_vme_associations.csv`
- `28_model_fit_diagnostics.csv`
- `28_models_excluded_numerical_instability.csv`

# Figure and GitHub-Linking System

Every analytical figure should have a deterministic number and descriptive filename.

Example:

```text
Figures/
├── Fig_01_Point_Composition_By_Dive.png
├── Fig_02_Substrate_Class_By_Dive.png
├── Fig_03_Richness_Depth.png
├── Fig_04_Shannon_Depth.png
├── Fig_05_NMDS_Latitude.png
├── Fig_06_NMDS_Substrate.png
├── Fig_07_CAP_Environment.png
└── Fig_08_VME_Depth.png
```

The README links directly to repository figures:

```markdown
![NMDS coloured by latitude](Figures/Fig_05_NMDS_Latitude.png)
```

The same section should link numerical outputs:

```markdown
[Download PERMANOVA table](Analyses/04_PERMANOVA/permanova_environment.csv)
```

## Recommended per-analysis GitHub pattern

```markdown
## PERMANOVA: environmental model

### Question
What environmental variables are associated with community composition?

### Input
- `Data/Final/point_annotations.csv`

### Script
- `Scripts/23_analysis_permanova.R`

### Model
`community ~ depth + temperature + substrate + relief`

### Result
[PERMANOVA result table](Analyses/04_PERMANOVA/permanova_environment.csv)

### Figure
![CAP environmental model](Figures/Fig_07_CAP_Environment.png)

<details>
<summary>How to interpret this result</summary>

Explain pseudo-F, R², p-values, dispersion test and limitations here.

</details>
```

This lets the README remain comprehensive without forcing every diagnostic table to be permanently expanded.

---

# Central Configuration

All user-facing options should be concentrated in `Scripts/00_config.R`, including:

- target point count;
- acceptable total-point range for inferential community models;
- minimum assigned biological points;
- biological-unit resolution;
- relative-abundance/square-root/other community transformation;
- Bray–Curtis distance;
- permutation number;
- candidate environmental variables;
- PERMANOVA formula terms;
- CAP terms;
- latitude palette;
- manual colour endpoints;
- point size/transparency;
- NMDS labels/ellipses;
- significance threshold;
- multiple-testing correction.

This means plot styling should be changed in **one obvious place**, not buried inside analysis code.

---

# Quality Control

At minimum, the formal QC module should produce:

## Input integrity

- presence of required files;
- required columns;
- duplicate/conflicting lookup keys;
- malformed rows.

## Join coverage

- `dive_id ↔ FullID` match rate;
- `filename ↔ FrameID` match rate;
- unmatched IDs written explicitly.

## Annotation design

- rows per frame;
- unique `annotation_id` values per frame;
- frames not meeting the expected ~100-point design;
- multiply exported annotation IDs;
- conflicting top-level labels;
- number of Biotic/Abiotic/Exclude/UNSURE points;
- excluded/uncertain image-quality proportion.

## Biological resolution

- proportion using `common_id_full`;
- proportion falling back to `common_id_mid`;
- proportion falling back to `common_id_short`;
- unresolved/multiple biological assignments.

## Environmental coverage

- frames with depth;
- frames with temperature;
- frames with coordinates;
- within-frame metadata consistency;
- ranges and impossible/outlying values.

## Statistical diagnostics

- sample size per group;
- predictor correlation/collinearity;
- model residual checks;
- NMDS stress;
- PERMDISP alongside categorical PERMANOVA;
- permutation restrictions;
- multiple-testing adjustment.

---

# Adding a New Analysis

A new analysis should be modular and should not edit upstream datasets.

Example: adding beta-diversity decomposition.

```text
Scripts/29_analysis_beta_diversity.R
Analyses/10_Beta_Diversity/
Figures/Fig_XX_Beta_Diversity_*.png
```

Then add a README section using the standard template:

```markdown
## Analysis title

### Ecological question
What are we asking?

### Sampling unit
Frame, dive, VME occurrence, etc.

### Response variable(s)
Exact column/matrix definition.

### Predictor(s)
Exact metadata fields.

### Data preparation
Filtering, transformations, standardisation.

### Statistical test
Name and implementation.

### Null hypothesis
State it explicitly.

### Assumptions / diagnostics
List checks.

### Output
Link tables/model diagnostics.

### Figure
Link numbered figure.

### Interpretation
Explain effect size/statistic/p-value and ecological meaning.

### Limitations
State pseudoreplication, missing metadata, observational inference, etc.
```

That template is the core mechanism that makes this repository expandable.

---

# Reproducibility and Interpretation

## Reproducibility principles

- source inputs remain unchanged;
- every transformation has a numbered script;
- every stage has a numbered output;
- all joins produce match statistics;
- random procedures use a documented random seed;
- package versions will ultimately be locked;
- figures are generated from scripts rather than edited manually;
- statistical outputs are saved as machine-readable CSV files where possible;
- model formulas and settings appear in output metadata/logs.

## Interpretation principles

> [!IMPORTANT]
> **Statistical significance is not ecological importance.** Report effect sizes, variance explained, fitted relationships and uncertainty alongside p-values.

> [!CAUTION]
> This is an observational image-analysis pipeline. Associations with depth, temperature, substrate, latitude or longitude do not by themselves demonstrate causation.

> [!CAUTION]
> Spatial and hierarchical dependence matters. Multiple frames within a dive/site may be more similar than frames from independent locations. Analyses must respect the sampling design through strata, aggregation, hierarchical models, or appropriately qualified inference.

> [!NOTE]
> Depth, temperature, latitude and other environmental predictors may themselves be correlated. Combined models should diagnose collinearity and avoid interpreting strongly confounded coefficients as independent causal effects.

---

# Development Status

| Component | Status |
|---|---|
| Canonical repository structure | Established |
| Canonical input filenames and prerequisite sheets | Established |
| Input validator and processing QC | Validated on worked example |
| Non-destructive Steps 01–07 processing | Validated on 55,096 source rows |
| Point/VME analytical split | Validated |
| Frame substrate characterisation | Validated |
| Shared R analysis helpers | Validated |
| Analysis 20 — point/frame QC | Validated |
| Analysis 21 — alpha diversity | Validated |
| Analysis 22 — Bray–Curtis/NMDS/envfit | Validated, including 3-D stress diagnostic |
| Analysis 23 — PERMANOVA/PERMDISP | Validated with dive-restricted permutations |
| Analysis 24 — CAP/dbRDA | Validated with VIF/alias diagnostics |
| Analysis 25 — environmental alpha-diversity models | Validated with dive hierarchy |
| Analysis 26 — hierarchical categorical comparisons | Validated; ineligible tests explicitly reported |
| Analysis 27 — focal-taxon responses | Validated with per-model convergence audit |
| Analysis 28 — VME occurrence | Validated with sampling-hierarchy and convergence audit |
| Compact worked-example result tables | Versioned in Git |
| GitHub-viewable figures | Fig. 01–12 versioned as PNG |
| Large reconstructable generated products | Ignored by Git and regenerated locally |
| Final clean full-pipeline rerun from canonical inputs | Pending release audit |
| Software lockfile / exact environment capture | Planned |
| Versioned release tag | Pending final clean rerun and documentation audit |

---

## Immediate Next Development Step

The analytical modules have been individually executed, reviewed and hardened on the worked example. The remaining release-stage tasks are:

1. run `./Scripts/run_full_pipeline.sh .` from the canonical prerequisite inputs;
2. confirm that the fresh full run reproduces the validated compact result tables and figures without new warnings/errors;
3. record the R/package environment used for the validated release;
4. perform a final README/output audit; and
5. create a versioned Git tag for the validated pipeline release.

No additional ecological hypothesis tests should be added merely because the current modules have completed. New analyses should be introduced only for a defined ecological question and should follow the eligibility/hierarchy rules documented above.
