# Scripts

This directory contains the executable stages of the BIIGLE benthic image-analysis pipeline.

## Working rule

Run scripts from the **repository root**, not from inside `Scripts/`.

For example:

```bash
cd ~/Documents/GitHub/BIIGLE-Benthic-Pipeline
./Scripts/run_processing_pipeline.sh .
```

## Central configuration

`00_config.R` contains user-editable analysis and plotting settings, including:

- expected point effort per frame;
- minimum total-point effort used in community models;
- minimum assigned biological points;
- rarefaction depth;
- community transformation;
- Bray-Curtis distance;
- NMDS dimensions and attempts;
- permutation count;
- environmental predictors;
- latitude colour mode and palette;
- figure dimensions;
- multiple-testing correction.

Edit this file before editing individual analysis scripts.

## Processing scripts

The processing pipeline is:

```text
00_validate_inputs.sh
        ↓
01_add_lat_long.sh
        ↓
02_add_top_level.sh
        ↓
03_add_taxonomy.sh
        ↓
04_add_frame_metadata.sh
        ↓
05_add_vme_taxonomy.sh
        ↓
06_add_frame_substrate_class.sh
        ↓
07_split_point_and_vme_data.sh
```

Run all stages with:

```bash
./Scripts/run_processing_pipeline.sh .
```

The scripts read the canonical prerequisite files from `Sheets/`, leave those files unchanged, and write numbered datasets to `Data/Intermediate/` and `Data/Final/`.

Steps 01–07 are non-destructive wrappers around versioned implementations retained in `internal/`. The `internal/` scripts are repository dependencies bundled with the pipeline; users normally run the numbered wrapper scripts rather than the implementation files directly.

## Analysis helpers

`01_analysis_helpers.R` contains shared functions for:

- resolving unique `filename + annotation_id` point observations;
- preventing multiply exported BIIGLE labels from becoming duplicate points;
- selecting the biological community unit using  
  `common_id_full → common_id_mid → common_id_short → label_name`;
- building frame summaries;
- building frame × biological-unit community matrices;
- calculating alpha-diversity metrics;
- applying configurable community transformations;
- saving numbered figures.

## Analysis scripts

```text
20_analysis_qc.R
21_analysis_alpha_diversity.R
22_analysis_nmds.R
23_analysis_permanova.R
24_analysis_cap_dbrda.R
25_analysis_environment.R
26_analysis_univariate.R
27_analysis_taxon_responses.R
28_analysis_vme.R
```

Run all analyses with:

```bash
./Scripts/run_analysis_pipeline.sh .
```

Or run one module, for example:

```bash
Rscript Scripts/22_analysis_nmds.R
```

## Full rerun

To rebuild all processed data and then execute every analysis:

```bash
./Scripts/run_full_pipeline.sh .
```

## R packages

Install required packages once with:

```bash
Rscript Scripts/install_R_packages.R
```

## Provenance

`internal/` contains versioned implementation scripts required by Steps 01–07. Do not run these directly for the canonical repository workflow unless debugging a processing stage.
