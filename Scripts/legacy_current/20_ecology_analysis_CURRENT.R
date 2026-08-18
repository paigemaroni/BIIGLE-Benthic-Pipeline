#!/usr/bin/env Rscript

# =============================================================================
# BIIGLE ECOLOGICAL ANALYSIS PIPELINE — POINT DATA + VME DATA
# =============================================================================
#
# This pipeline is designed for the two-table BIIGLE workflow:
#
#   1. point_annotations.csv
#      Random-point annotation data. Frames should contain approximately
#      TARGET_POINTS_PER_FRAME (normally 100) unique annotation points.
#      This file drives cover, substrate, diversity, community, NMDS, GAM and
#      PERMANOVA analyses.
#
#   2. vme_annotations.csv
#      VME observations (typically polygon annotations) recorded only when a
#      VME was observed in a frame. These records are analysed separately as
#      VME occurrence/presence evidence and NEVER enter the 100-point
#      denominator, biological relative abundance, richness or Bray-Curtis
#      community matrices.
#
# CURRENT CANONICAL BIOLOGICAL FIELDS
# -----------------------------------
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
# CURRENT CANONICAL ABIOTIC / HABITAT FIELDS
# ------------------------------------------
#   releif       (spelling preserved from TAXON_LIST.csv)
#   substrate
#   type
#   size
#
# CURRENT ENVIRONMENTAL FIELDS
# ----------------------------
#   video_time
#   depth_m
#   temperature
#
# POINT-DATA COMMUNITY UNIT
# -------------------------
#   common_id_full -> common_id_mid -> common_id_short -> label_name
#
# IMPORTANT INTERPRETATION
# ------------------------
# - Point percentages and diversity come ONLY from point_annotations.csv.
# - Exclude and UNSURE remain QC categories and are excluded from valid
#   ecological denominators.
# - Abiotic random points contribute habitat/substrate information but never
#   enter the biological community matrix.
# - VME polygons are not random points and are therefore not treated as point
#   abundance. VME data are summarised as annotation occurrence, frame-level
#   presence/absence, VME-type richness/occurrence, and dive-level occurrence.
# - A frame present in point_annotations.csv but absent from vme_annotations.csv
#   is treated as having no recorded VME observation for this dataset.
# - VME frames that cannot be matched to a point-data frame are retained in
#   VME-specific outputs and reported for review.
#
# QUESTIONS
# ---------
# Q1. Does image quality (Exclude percentage) vary with depth or temperature?
# Q2. Does biotic cover vary with depth, temperature, dive, or substrate?
# Q3. Does standardised biological-group richness/diversity vary with depth
#     or temperature?
# Q4. Does point-community composition turn over among frames and dives?
# Q5. Are depth, temperature, and substrate associated with point-community
#     composition?
# Q6. Which biological groups have the strongest depth or temperature
#     associations in the point data?
# Q7. Which frames and dives are candidate biodiversity hotspots from the
#     random-point data?
# Q8. Where are VMEs recorded, which VME types occur, and how does VME presence
#     overlap with point-data frames and environmental gradients?
#
# RUN
# ---
# With the default filenames in the working directory:
#
#   Rscript biigle_ecology_pipeline.R
#
# Or explicitly:
#
#   Rscript biigle_ecology_pipeline.R \
#       point_annotations.csv \
#       vme_annotations.csv \
#       biigle_ecology_analysis
#
# =============================================================================

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

POINT_INPUT_DEFAULT <- "point_annotations.csv"
VME_INPUT_DEFAULT <- "vme_annotations.csv"
OUTPUT_DEFAULT <- "biigle_ecology_analysis"

# Required BIIGLE identity columns
DIVE_ID_COL <- "dive_id"
FRAME_COL <- "filename"
POINT_ID_COL <- "annotation_id"
CLASS_COL <- "top_level"
HIERARCHY_COL <- "label_hierarchy"
LABEL_NAME_COL <- "label_name"

# Standardised TAXON_LIST biological fields
CPC_CODES_COL <- "cpc_codes"
KINGDOM_COL <- "kingdom"
PHYLUM_COL <- "phylum"
TAXON_CLASS_COL <- "class"
ORDER_COL <- "order"
FAMILY_COL <- "family"
TAXONOMIC_RESOLUTION_COL <- "taxonomic_resolution"
MORPHOLOGY_COL <- "morphology"
COMMON_ID_SHORT_COL <- "common_id_short"
COMMON_ID_MID_COL <- "common_id_mid"
COMMON_ID_FULL_COL <- "common_id_full"

# Standardised TAXON_LIST abiotic/habitat fields.
# NOTE: "releif" intentionally matches the current CSV header.
RELIEF_COL <- "releif"
SUBSTRATE_COL <- "substrate"
SUBSTRATE_TYPE_COL <- "type"
SUBSTRATE_SIZE_COL <- "size"

# Frame-level environmental metadata
VIDEO_TIME_COL <- "video_time"
DEPTH_COL <- "depth_m"
TEMPERATURE_COL <- "temperature"
LATITUDE_COL <- "lat"
LONGITUDE_COL <- "long"

# Optional legacy time column.
TIME_COL <- "time"

# Expected unique random points per frame in point_annotations.csv.
TARGET_POINTS_PER_FRAME <- 100

# Minimum assigned biological random points for frame-level diversity/ordination.
MIN_BIOTIC_POINTS <- 10

# Standardised richness is rarefied to this many biological random points.
RAREFACTION_N <- 10

# Number of common point-data groups shown in response plots and heatmaps.
TOP_N_GROUPS <- 12

# Exploratory models are generated when data coverage permits.
RUN_EXPLORATORY_GAMS <- TRUE
RUN_PERMANOVA <- TRUE

N_PERMUTATIONS <- 999
RANDOM_SEED <- 20260810

set.seed(RANDOM_SEED)

# =============================================================================
# 2. PACKAGES
# =============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "vegan",
  "mgcv"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them with:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(vegan)
library(mgcv)


# =============================================================================
# 3. INPUT AND OUTPUT PATHS
# =============================================================================

arguments <- commandArgs(trailingOnly = TRUE)

point_input_csv <- if (length(arguments) >= 1) {
  arguments[[1]]
} else {
  POINT_INPUT_DEFAULT
}

vme_input_csv <- if (length(arguments) >= 2) {
  arguments[[2]]
} else {
  VME_INPUT_DEFAULT
}

output_root <- if (length(arguments) >= 3) {
  arguments[[3]]
} else {
  OUTPUT_DEFAULT
}

if (!file.exists(point_input_csv)) {
  stop(
    paste0("Point-data input file not found: ", point_input_csv),
    call. = FALSE
  )
}

if (!file.exists(vme_input_csv)) {
  stop(
    paste0("VME input file not found: ", vme_input_csv),
    call. = FALSE
  )
}

directories <- list(
  qc = file.path(output_root, "01_quality_control"),
  composition = file.path(output_root, "02_frame_composition"),
  community = file.path(output_root, "03_community_matrices"),
  diversity = file.path(output_root, "04_diversity"),
  environment = file.path(output_root, "05_environment"),
  ordination = file.path(output_root, "06_ordination"),
  groups = file.path(output_root, "07_group_responses"),
  vme = file.path(output_root, "08_vme_occurrence"),
  models = file.path(output_root, "09_models"),
  summaries = file.path(output_root, "10_summaries"),
  logs = file.path(output_root, "11_logs")
)

invisible(
  lapply(
    directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

assert_columns <- function(data, columns, context = "input data") {
  missing <- setdiff(columns, names(data))

  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required columns in ",
        context,
        ": ",
        paste(missing, collapse = ", "),
        "\n\nColumns found:\n",
        paste(names(data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
}


column_or_na <- function(data, column_name, numeric = FALSE) {
  if (column_name %in% names(data)) {
    return(data[[column_name]])
  }

  if (numeric) {
    return(rep(NA_real_, nrow(data)))
  }

  rep(NA_character_, nrow(data))
}


clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}


safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}


median_or_na <- function(x) {
  x <- safe_numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(x, na.rm = TRUE)
}


range_or_na <- function(x) {
  x <- safe_numeric(x)

  if (all(is.na(x))) {
    return(NA_real_)
  }

  diff(range(x, na.rm = TRUE))
}


first_non_missing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != ""]

  if (length(x) == 0) {
    return(NA_character_)
  }

  as.character(x[[1]])
}


collapse_unique <- function(x) {
  values <- sort(unique(clean_character(x)))
  values <- values[!is.na(values)]

  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(values, collapse = " | ")
}


count_unique <- function(x) {
  values <- unique(clean_character(x))
  sum(!is.na(values))
}


weighted_mean_or_na <- function(x, weights) {
  x <- safe_numeric(x)
  weights <- safe_numeric(weights)

  keep <- !is.na(x) & !is.na(weights) & weights > 0

  if (!any(keep)) {
    return(NA_real_)
  }

  weighted.mean(x[keep], weights[keep])
}


scale_safe <- function(x) {
  x <- safe_numeric(x)

  if (
    all(is.na(x)) ||
      is.na(sd(x, na.rm = TRUE)) ||
      sd(x, na.rm = TRUE) == 0
  ) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
}


standardise_top_level <- function(x) {
  x_clean <- str_to_upper(str_trim(as.character(x)))

  case_when(
    str_detect(x_clean, "^BIOTIC$") ~ "Biotic",
    str_detect(x_clean, "^ABIOTIC$") ~ "Abiotic",
    str_detect(x_clean, "^EXCLUDE$") ~ "Exclude",
    str_detect(x_clean, "^UNSURE$") ~ "Unsure",
    str_detect(x_clean, "^VME$") ~ "VME",
    is.na(x_clean) | x_clean == "" ~ "Unclassified",
    TRUE ~ "Other"
  )
}


split_hierarchy_vector <- function(x) {
  str_split(
    as.character(x),
    "\\s*>\\s*"
  )
}


hierarchy_entry <- function(parts_list, position) {
  vapply(
    parts_list,
    function(parts) {
      parts <- trimws(parts)
      parts <- parts[parts != ""]

      if (length(parts) >= position) {
        return(parts[[position]])
      }

      NA_character_
    },
    FUN.VALUE = character(1)
  )
}


value_after_token <- function(parts_list, token) {
  token_clean <- str_to_lower(token)

  vapply(
    parts_list,
    function(parts) {
      parts <- trimws(parts)
      parts <- parts[parts != ""]

      index <- which(str_to_lower(parts) == token_clean)

      if (
        length(index) > 0 &&
          index[[1]] < length(parts)
      ) {
        return(parts[[index[[1]] + 1]])
      }

      NA_character_
    },
    FUN.VALUE = character(1)
  )
}


first_matching_value <- function(parts_list, accepted_values) {
  accepted_clean <- str_to_lower(accepted_values)

  vapply(
    parts_list,
    function(parts) {
      parts <- trimws(parts)
      parts <- parts[parts != ""]

      match_index <- which(
        str_to_lower(parts) %in% accepted_clean
      )

      if (length(match_index) > 0) {
        return(parts[[match_index[[1]]]])
      }

      NA_character_
    },
    FUN.VALUE = character(1)
  )
}


save_plot <- function(
    plot,
    filename_stub,
    directory,
    width = 10,
    height = 7,
    dpi = 300
) {
  ggsave(
    filename = file.path(
      directory,
      paste0(filename_stub, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(
      directory,
      paste0(filename_stub, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}


write_text <- function(text, path) {
  writeLines(
    text,
    con = path,
    useBytes = TRUE
  )
}


make_environment_plot <- function(
    data,
    predictor,
    response,
    title,
    y_label,
    output_stub,
    output_directory,
    percent_axis = FALSE
) {
  plot_data <- data %>%
    filter(
      !is.na(.data[[predictor]]),
      !is.na(.data[[response]])
    )

  if (
    nrow(plot_data) < 4 ||
      n_distinct(plot_data[[predictor]]) < 3
  ) {
    return(invisible(NULL))
  }

  plot_object <- ggplot(
    plot_data,
    aes(
      x = .data[[predictor]],
      y = .data[[response]],
      colour = dive_id_analysis
    )
  ) +
    geom_point(
      alpha = 0.75,
      size = 2
    ) +
    labs(
      title = title,
      x = ifelse(
        predictor == "depth_m",
        "Depth (m; positive downward)",
        "Temperature"
      ),
      y = y_label,
      colour = "Dive ID"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  if (
    nrow(plot_data) >= 12 &&
      n_distinct(plot_data[[predictor]]) >= 6
  ) {
    plot_object <- plot_object +
      geom_smooth(
        data = plot_data,
        mapping = aes(
          x = .data[[predictor]],
          y = .data[[response]],
          group = 1
        ),
        inherit.aes = FALSE,
        method = "gam",
        formula = y ~ s(x, bs = "cs", k = 4),
        se = TRUE,
        colour = "black",
        linewidth = 0.8
      )
  }

  if (percent_axis) {
    plot_object <- plot_object +
      scale_y_continuous(
        limits = c(0, 100),
        expand = expansion(mult = c(0.02, 0.05))
      )
  }

  save_plot(
    plot_object,
    output_stub,
    output_directory,
    width = 9,
    height = 7
  )

  invisible(plot_object)
}


safe_spearman <- function(x, y) {
  keep <- complete.cases(x, y)

  if (
    sum(keep) < 6 ||
      n_distinct(x[keep]) < 3 ||
      n_distinct(y[keep]) < 2
  ) {
    return(
      data.frame(
        rho = NA_real_,
        p_value = NA_real_,
        n = sum(keep)
      )
    )
  }

  result <- suppressWarnings(
    cor.test(
      x[keep],
      y[keep],
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(
    rho = unname(result$estimate),
    p_value = result$p.value,
    n = sum(keep)
  )
}


# =============================================================================
# 5. READ POINT + VME DATA AND IDENTIFY COMMUNITY-UNIT MODE
# =============================================================================

point_raw_data <- read_csv(
  point_input_csv,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "N/A", "NULL")
)

vme_raw_data <- read_csv(
  vme_input_csv,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "N/A", "NULL")
)

assert_columns(
  point_raw_data,
  c(
    DIVE_ID_COL,
    FRAME_COL,
    POINT_ID_COL,
    CLASS_COL,
    HIERARCHY_COL,
    COMMON_ID_SHORT_COL,
    COMMON_ID_MID_COL,
    COMMON_ID_FULL_COL,
    RELIEF_COL,
    SUBSTRATE_COL,
    SUBSTRATE_TYPE_COL,
    SUBSTRATE_SIZE_COL,
    VIDEO_TIME_COL,
    DEPTH_COL,
    TEMPERATURE_COL
  ),
  context = "point_annotations.csv"
)

assert_columns(
  vme_raw_data,
  c(
    DIVE_ID_COL,
    FRAME_COL,
    POINT_ID_COL,
    CLASS_COL,
    HIERARCHY_COL,
    CPC_CODES_COL,
    KINGDOM_COL,
    PHYLUM_COL,
    TAXON_CLASS_COL,
    ORDER_COL,
    FAMILY_COL,
    TAXONOMIC_RESOLUTION_COL,
    MORPHOLOGY_COL,
    COMMON_ID_SHORT_COL,
    COMMON_ID_MID_COL,
    COMMON_ID_FULL_COL
  ),
  context = "vme_annotations.csv"
)

# Keep raw_data as a point-data alias so the existing point-analysis branch
# remains easy to audit. VME data are handled separately below.
raw_data <- point_raw_data

hierarchy_raw <- clean_character(
  column_or_na(raw_data, HIERARCHY_COL)
)

hierarchy_parts <- split_hierarchy_vector(hierarchy_raw)

class_original <- standardise_top_level(
  column_or_na(raw_data, CLASS_COL)
)

# The point file must not contain VME polygons. Stopping here prevents an
# accidental double-count if the split step was not completed correctly.
point_vme_rows <- sum(class_original == "VME", na.rm = TRUE)

if (point_vme_rows > 0) {
  stop(
    paste0(
      "point_annotations.csv still contains ",
      point_vme_rows,
      " VME row(s). Split VME annotations out before running this pipeline."
    ),
    call. = FALSE
  )
}

# The VME file is expected to contain only VME observations.
vme_class_check <- standardise_top_level(
  column_or_na(vme_raw_data, CLASS_COL)
)

non_vme_in_vme_file <- sum(
  is.na(vme_class_check) | vme_class_check != "VME"
)

if (non_vme_in_vme_file > 0) {
  stop(
    paste0(
      "vme_annotations.csv contains ",
      non_vme_in_vme_file,
      " non-VME row(s). The VME file should contain top_level == VME only."
    ),
    call. = FALSE
  )
}

# Standardised common identifiers for point-community analysis.
common_id_short_raw <- clean_character(
  column_or_na(raw_data, COMMON_ID_SHORT_COL)
)

common_id_mid_raw <- clean_character(
  column_or_na(raw_data, COMMON_ID_MID_COL)
)

common_id_full_raw <- clean_character(
  column_or_na(raw_data, COMMON_ID_FULL_COL)
)

label_name_raw <- clean_character(
  column_or_na(raw_data, LABEL_NAME_COL)
)

community_unit_mode <- paste0(
  "point-data standardised common ID: ",
  COMMON_ID_FULL_COL,
  " -> ",
  COMMON_ID_MID_COL,
  " -> ",
  COMMON_ID_SHORT_COL,
  " -> label_name fallback"
)

# =============================================================================
# 6. CREATE ROW-LEVEL ANALYSIS FIELDS
# =============================================================================

# ---------------------------------------------------------------------------
# Biological classification from the current TAXON_LIST fields
# ---------------------------------------------------------------------------

cpc_codes_raw <- clean_character(
  column_or_na(raw_data, CPC_CODES_COL)
)

kingdom_raw <- clean_character(
  column_or_na(raw_data, KINGDOM_COL)
)

phylum_raw <- clean_character(
  column_or_na(raw_data, PHYLUM_COL)
)

taxon_class_raw <- clean_character(
  column_or_na(raw_data, TAXON_CLASS_COL)
)

order_raw <- clean_character(
  column_or_na(raw_data, ORDER_COL)
)

family_raw <- clean_character(
  column_or_na(raw_data, FAMILY_COL)
)

taxonomic_resolution_raw <- clean_character(
  column_or_na(raw_data, TAXONOMIC_RESOLUTION_COL)
)

morphology_raw <- clean_character(
  column_or_na(raw_data, MORPHOLOGY_COL)
)

# Standardised biological hierarchy used for ecological analyses.
# Restrict these fields to Biotic random-point rows so abiotic common IDs can never
# enter biological richness or community calculations. VME data are separate.
broad_group_row <- case_when(
  class_original == "Biotic" ~
    common_id_short_raw,
  TRUE ~ NA_character_
)

mid_group_row <- case_when(
  class_original == "Biotic" ~
    coalesce(
      common_id_mid_raw,
      common_id_short_raw
    ),
  TRUE ~ NA_character_
)

fine_group_row <- case_when(
  class_original == "Biotic" ~
    coalesce(
      common_id_full_raw,
      common_id_mid_raw,
      common_id_short_raw
    ),
  TRUE ~ NA_character_
)

community_unit_row <- case_when(
  class_original == "Biotic" ~
    coalesce(
      common_id_full_raw,
      common_id_mid_raw,
      common_id_short_raw,
      label_name_raw
    ),
  TRUE ~ NA_character_
)

community_unit_source_row <- case_when(
  class_original == "Biotic" &
    !is.na(common_id_full_raw) ~ "common_id_full",
  class_original == "Biotic" &
    !is.na(common_id_mid_raw) ~ "common_id_mid",
  class_original == "Biotic" &
    !is.na(common_id_short_raw) ~ "common_id_short",
  class_original == "Biotic" &
    !is.na(label_name_raw) ~ "label_name_fallback",
  TRUE ~ NA_character_
)


# ---------------------------------------------------------------------------
# Abiotic / habitat classification from the current TAXON_LIST fields
# ---------------------------------------------------------------------------

relief_raw <- clean_character(
  column_or_na(raw_data, RELIEF_COL)
)

substrate_raw <- clean_character(
  column_or_na(raw_data, SUBSTRATE_COL)
)

substrate_method_raw <- clean_character(
  column_or_na(raw_data, SUBSTRATE_TYPE_COL)
)

substrate_size_raw <- clean_character(
  column_or_na(raw_data, SUBSTRATE_SIZE_COL)
)

# Retain hierarchy fallbacks only for older/partially enriched rows.
substrate_from_hierarchy <- first_matching_value(
  hierarchy_parts,
  c("Hard", "Soft")
)

substrate_size_from_hierarchy <- value_after_token(
  hierarchy_parts,
  "Wentworth scale"
)

relief_from_hierarchy <- value_after_token(
  hierarchy_parts,
  "Habitat Relief"
)

relief_row <- case_when(
  class_original == "Abiotic" ~
    coalesce(
      relief_raw,
      relief_from_hierarchy
    ),
  TRUE ~ NA_character_
)

substrate_type_row <- case_when(
  class_original == "Abiotic" ~
    coalesce(
      substrate_raw,
      substrate_from_hierarchy
    ),
  TRUE ~ NA_character_
)

substrate_method_row <- case_when(
  class_original == "Abiotic" ~
    substrate_method_raw,
  TRUE ~ NA_character_
)

substrate_size_row <- case_when(
  class_original == "Abiotic" ~
    coalesce(
      substrate_size_raw,
      substrate_size_from_hierarchy
    ),
  TRUE ~ NA_character_
)

# "substrate_detail" is kept for compatibility with the existing downstream
# plots. It uses the finest available abiotic descriptor.
substrate_detail_row <- case_when(
  class_original == "Abiotic" ~
    coalesce(
      substrate_size_row,
      substrate_method_row,
      substrate_type_row
    ),
  TRUE ~ NA_character_
)


class_ecological <- case_when(
  class_original == "Unsure" ~ "Unsure",
  TRUE ~ class_original
)

row_data <- raw_data %>%
  mutate(
    dive_id_analysis = clean_character(
      .data[[DIVE_ID_COL]]
    ),
    frame_id_analysis = clean_character(
      .data[[FRAME_COL]]
    ),
    point_id_analysis = clean_character(
      .data[[POINT_ID_COL]]
    ),
    sample_id_analysis = paste(
      dive_id_analysis,
      frame_id_analysis,
      sep = "__"
    ),
    hierarchy_analysis = hierarchy_raw,
    class_original_analysis = class_original,
    class_ecological_analysis = class_ecological,
    cpc_codes_row = cpc_codes_raw,
    kingdom_row = kingdom_raw,
    phylum_row = phylum_raw,
    taxon_class_row = taxon_class_raw,
    order_row = order_raw,
    family_row = family_raw,
    taxonomic_resolution_row = taxonomic_resolution_raw,
    morphology_row = morphology_raw,
    common_id_short_row = common_id_short_raw,
    common_id_mid_row = common_id_mid_raw,
    common_id_full_row = common_id_full_raw,
    broad_group_row = broad_group_row,
    mid_group_row = mid_group_row,
    fine_group_row = fine_group_row,
    community_unit_row = community_unit_row,
    community_unit_source_row = community_unit_source_row,
    relief_row = relief_row,
    substrate_type_row = substrate_type_row,
    substrate_method_row = substrate_method_row,
    substrate_size_row = substrate_size_row,
    substrate_detail_row = substrate_detail_row,
    depth_raw_m_analysis = safe_numeric(
      column_or_na(raw_data, DEPTH_COL, numeric = TRUE)
    ),
    # Uploaded depth values are negative below the surface. Use positive
    # magnitude for ecological interpretation and plotting.
    depth_m_analysis = abs(depth_raw_m_analysis),
    temperature_analysis = safe_numeric(
      column_or_na(
        raw_data,
        TEMPERATURE_COL,
        numeric = TRUE
      )
    ),
    latitude_analysis = safe_numeric(
      column_or_na(
        raw_data,
        LATITUDE_COL,
        numeric = TRUE
      )
    ),
    longitude_analysis = safe_numeric(
      column_or_na(
        raw_data,
        LONGITUDE_COL,
        numeric = TRUE
      )
    ),
    time_analysis = clean_character(
      column_or_na(raw_data, TIME_COL)
    ),
    video_time_analysis = clean_character(
      column_or_na(raw_data, VIDEO_TIME_COL)
    )
  )


# =============================================================================
# 7. MISSING-KEY AND RAW-CLASS AUDITS
# =============================================================================

missing_key_rows <- row_data %>%
  filter(
    is.na(dive_id_analysis) |
      is.na(frame_id_analysis) |
      is.na(point_id_analysis)
  )

write_csv(
  missing_key_rows,
  file.path(
    directories$qc,
    "rows_missing_dive_frame_or_point_id.csv"
  )
)

row_data_valid <- row_data %>%
  filter(
    !is.na(dive_id_analysis),
    !is.na(frame_id_analysis),
    !is.na(point_id_analysis)
  )

class_audit <- row_data_valid %>%
  count(
    class_original_analysis,
    class_ecological_analysis,
    name = "n_rows",
    sort = TRUE
  )

write_csv(
  class_audit,
  file.path(
    directories$qc,
    "top_level_classification_audit.csv"
  )
)


# =============================================================================
# 8. RESOLVE MULTIPLY LABELLED ANNOTATION POINTS
# =============================================================================
#
# A single BIIGLE annotation_id can contain more than one exported label row.
# The pipeline resolves these at the unique point level.
#
# - Repeated rows with one ecological class remain one point.
# - Points assigned to multiple ecological classes are marked Conflict.
# - Biological points with multiple different standardised common-ID units are retained
#   for cover but excluded from diversity/community matrices.
# =============================================================================

primary_classes <- c(
  "Biotic",
  "Abiotic",
  "Exclude",
  "Unsure"
)

point_data <- row_data_valid %>%
  group_by(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    point_id_analysis
  ) %>%
  summarise(
    exported_rows_per_point = n(),
    primary_class_set = collapse_unique(
      class_ecological_analysis[
        class_ecological_analysis %in% primary_classes
      ]
    ),
    n_primary_classes = count_unique(
      class_ecological_analysis[
        class_ecological_analysis %in% primary_classes
      ]
    ),
    cpc_codes_set = collapse_unique(
      cpc_codes_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    kingdom_set = collapse_unique(
      kingdom_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    phylum_set = collapse_unique(
      phylum_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    taxon_class_set = collapse_unique(
      taxon_class_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    order_set = collapse_unique(
      order_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    family_set = collapse_unique(
      family_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    taxonomic_resolution_set = collapse_unique(
      taxonomic_resolution_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    morphology_set = collapse_unique(
      morphology_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    common_id_short_set = collapse_unique(
      common_id_short_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    common_id_mid_set = collapse_unique(
      common_id_mid_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    common_id_full_set = collapse_unique(
      common_id_full_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    broad_group_set = collapse_unique(
      broad_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    n_broad_groups = count_unique(
      broad_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    mid_group_set = collapse_unique(
      mid_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    n_mid_groups = count_unique(
      mid_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    fine_group_set = collapse_unique(
      fine_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    n_fine_groups = count_unique(
      fine_group_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    community_unit_set = collapse_unique(
      community_unit_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    n_community_units = count_unique(
      community_unit_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    community_unit_source_set = collapse_unique(
      community_unit_source_row[
        class_ecological_analysis == "Biotic"
      ]
    ),
    substrate_type_set = collapse_unique(
      substrate_type_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    n_substrate_types = count_unique(
      substrate_type_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    substrate_detail_set = collapse_unique(
      substrate_detail_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    n_substrate_details = count_unique(
      substrate_detail_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    substrate_method_set = collapse_unique(
      substrate_method_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    n_substrate_methods = count_unique(
      substrate_method_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    substrate_size_set = collapse_unique(
      substrate_size_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    n_substrate_sizes = count_unique(
      substrate_size_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    relief_set = collapse_unique(
      relief_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    n_relief_values = count_unique(
      relief_row[
        class_ecological_analysis == "Abiotic"
      ]
    ),
    depth_raw_m = median_or_na(
      depth_raw_m_analysis
    ),
    depth_m = median_or_na(
      depth_m_analysis
    ),
    temperature = median_or_na(
      temperature_analysis
    ),
    latitude = median_or_na(
      latitude_analysis
    ),
    longitude = median_or_na(
      longitude_analysis
    ),
    time = first_non_missing(
      time_analysis
    ),
    video_time = first_non_missing(
      video_time_analysis
    ),
    .groups = "drop"
  ) %>%
  mutate(
    point_class = case_when(
      n_primary_classes == 1 ~ primary_class_set,
      n_primary_classes == 0 ~ "Other",
      TRUE ~ "Conflict"
    ),
    cpc_codes = if_else(
      point_class == "Biotic",
      cpc_codes_set,
      NA_character_
    ),
    kingdom = if_else(
      point_class == "Biotic",
      kingdom_set,
      NA_character_
    ),
    phylum = if_else(
      point_class == "Biotic",
      phylum_set,
      NA_character_
    ),
    taxon_class = if_else(
      point_class == "Biotic",
      taxon_class_set,
      NA_character_
    ),
    order = if_else(
      point_class == "Biotic",
      order_set,
      NA_character_
    ),
    family = if_else(
      point_class == "Biotic",
      family_set,
      NA_character_
    ),
    taxonomic_resolution = if_else(
      point_class == "Biotic",
      taxonomic_resolution_set,
      NA_character_
    ),
    morphology = if_else(
      point_class == "Biotic",
      morphology_set,
      NA_character_
    ),
    common_id_short = if_else(
      point_class == "Biotic",
      common_id_short_set,
      NA_character_
    ),
    common_id_mid = if_else(
      point_class == "Biotic",
      common_id_mid_set,
      NA_character_
    ),
    common_id_full = if_else(
      point_class == "Biotic",
      common_id_full_set,
      NA_character_
    ),
    broad_group = if_else(
      point_class == "Biotic" &
        n_broad_groups == 1,
      broad_group_set,
      NA_character_
    ),
    mid_group = if_else(
      point_class == "Biotic" &
        n_mid_groups == 1,
      mid_group_set,
      NA_character_
    ),
    fine_group = if_else(
      point_class == "Biotic" &
        n_fine_groups == 1,
      fine_group_set,
      NA_character_
    ),
    community_unit = if_else(
      point_class == "Biotic" &
        n_community_units == 1,
      community_unit_set,
      NA_character_
    ),
    multiple_community_units =
      point_class == "Biotic" &
      n_community_units > 1,
    substrate_type = if_else(
      point_class == "Abiotic" &
        n_substrate_types == 1,
      substrate_type_set,
      NA_character_
    ),
    substrate_detail = if_else(
      point_class == "Abiotic" &
        n_substrate_details == 1,
      substrate_detail_set,
      NA_character_
    ),
    substrate_method = if_else(
      point_class == "Abiotic" &
        n_substrate_methods == 1,
      substrate_method_set,
      NA_character_
    ),
    substrate_size = if_else(
      point_class == "Abiotic" &
        n_substrate_sizes == 1,
      substrate_size_set,
      NA_character_
    ),
    relief = if_else(
      point_class == "Abiotic" &
        n_relief_values == 1,
      relief_set,
      NA_character_
    )
  )

write_csv(
  point_data,
  file.path(
    directories$qc,
    "resolved_unique_points.csv"
  )
)

conflicting_points <- point_data %>%
  filter(point_class == "Conflict")

write_csv(
  conflicting_points,
  file.path(
    directories$qc,
    "points_with_conflicting_top_levels.csv"
  )
)

multiple_unit_points <- point_data %>%
  filter(multiple_community_units)

write_csv(
  multiple_unit_points,
  file.path(
    directories$qc,
    "biotic_points_with_multiple_community_units.csv"
  )
)

duplicate_export_points <- point_data %>%
  filter(exported_rows_per_point > 1)

write_csv(
  duplicate_export_points,
  file.path(
    directories$qc,
    "multiply_labelled_annotation_points.csv"
  )
)


# =============================================================================
# 8B. VME OCCURRENCE DATA — SEPARATE FROM RANDOM-POINT ECOLOGY
# =============================================================================
#
# VME records are typically polygons/targeted observations. They are therefore
# summarised as occurrence/presence evidence only and never enter the random-
# point denominator or point-community matrices.
# =============================================================================

vme_common_id_short_raw <- clean_character(
  column_or_na(vme_raw_data, COMMON_ID_SHORT_COL)
)

vme_common_id_mid_raw <- clean_character(
  column_or_na(vme_raw_data, COMMON_ID_MID_COL)
)

vme_common_id_full_raw <- clean_character(
  column_or_na(vme_raw_data, COMMON_ID_FULL_COL)
)

vme_label_name_raw <- clean_character(
  column_or_na(vme_raw_data, LABEL_NAME_COL)
)

vme_unit_row <- coalesce(
  vme_common_id_full_raw,
  vme_common_id_mid_raw,
  vme_common_id_short_raw,
  vme_label_name_raw
)

vme_data <- vme_raw_data %>%
  mutate(
    dive_id_analysis = clean_character(
      .data[[DIVE_ID_COL]]
    ),
    frame_id_analysis = clean_character(
      .data[[FRAME_COL]]
    ),
    point_id_analysis = clean_character(
      .data[[POINT_ID_COL]]
    ),
    sample_id_analysis = paste(
      dive_id_analysis,
      frame_id_analysis,
      sep = "__"
    ),
    vme_unit = vme_unit_row,
    cpc_codes_vme = clean_character(
      column_or_na(vme_raw_data, CPC_CODES_COL)
    ),
    kingdom_vme = clean_character(
      column_or_na(vme_raw_data, KINGDOM_COL)
    ),
    phylum_vme = clean_character(
      column_or_na(vme_raw_data, PHYLUM_COL)
    ),
    class_vme = clean_character(
      column_or_na(vme_raw_data, TAXON_CLASS_COL)
    ),
    order_vme = clean_character(
      column_or_na(vme_raw_data, ORDER_COL)
    ),
    family_vme = clean_character(
      column_or_na(vme_raw_data, FAMILY_COL)
    ),
    taxonomic_resolution_vme = clean_character(
      column_or_na(vme_raw_data, TAXONOMIC_RESOLUTION_COL)
    ),
    morphology_vme = clean_character(
      column_or_na(vme_raw_data, MORPHOLOGY_COL)
    ),
    common_id_short_vme = vme_common_id_short_raw,
    common_id_mid_vme = vme_common_id_mid_raw,
    common_id_full_vme = vme_common_id_full_raw,
    depth_raw_m_vme = safe_numeric(
      column_or_na(vme_raw_data, DEPTH_COL, numeric = TRUE)
    ),
    depth_m_vme = abs(depth_raw_m_vme),
    temperature_vme = safe_numeric(
      column_or_na(vme_raw_data, TEMPERATURE_COL, numeric = TRUE)
    ),
    latitude_vme = safe_numeric(
      column_or_na(vme_raw_data, LATITUDE_COL, numeric = TRUE)
    ),
    longitude_vme = safe_numeric(
      column_or_na(vme_raw_data, LONGITUDE_COL, numeric = TRUE)
    ),
    video_time_vme = clean_character(
      column_or_na(vme_raw_data, VIDEO_TIME_COL)
    )
  )

vme_missing_key_rows <- vme_data %>%
  filter(
    is.na(dive_id_analysis) |
      is.na(frame_id_analysis) |
      is.na(point_id_analysis)
  )

write_csv(
  vme_missing_key_rows,
  file.path(
    directories$qc,
    "vme_rows_missing_dive_frame_or_annotation_id.csv"
  )
)

vme_data_valid <- vme_data %>%
  filter(
    !is.na(dive_id_analysis),
    !is.na(frame_id_analysis),
    !is.na(point_id_analysis)
  )

# One VME annotation is one occurrence record, even if an exported annotation
# contains more than one label row. Taxon/type summaries retain the distinct
# VME units attached to that annotation.
vme_annotations_unique <- vme_data_valid %>%
  group_by(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    point_id_analysis
  ) %>%
  summarise(
    exported_rows_per_vme_annotation = n(),
    vme_unit_set = collapse_unique(vme_unit),
    n_vme_units = count_unique(vme_unit),
    depth_m_vme = median_or_na(depth_m_vme),
    temperature_vme = median_or_na(temperature_vme),
    latitude_vme = median_or_na(latitude_vme),
    longitude_vme = median_or_na(longitude_vme),
    video_time_vme = first_non_missing(video_time_vme),
    .groups = "drop"
  )

write_csv(
  vme_annotations_unique,
  file.path(
    directories$vme,
    "vme_unique_annotations.csv"
  )
)

# Frame x VME-type occurrence table. A VME type is counted once per annotation.
vme_frame_type_occurrence <- vme_data_valid %>%
  filter(!is.na(vme_unit)) %>%
  distinct(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    point_id_analysis,
    vme_unit
  ) %>%
  count(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    vme_unit,
    name = "n_vme_annotations"
  ) %>%
  mutate(presence = 1L)

write_csv(
  vme_frame_type_occurrence,
  file.path(
    directories$vme,
    "vme_frame_type_occurrence.csv"
  )
)

vme_taxonomy <- vme_data_valid %>%
  filter(!is.na(vme_unit)) %>%
  group_by(vme_unit) %>%
  summarise(
    cpc_codes = collapse_unique(cpc_codes_vme),
    kingdom = collapse_unique(kingdom_vme),
    phylum = collapse_unique(phylum_vme),
    class = collapse_unique(class_vme),
    order = collapse_unique(order_vme),
    family = collapse_unique(family_vme),
    taxonomic_resolution = collapse_unique(taxonomic_resolution_vme),
    morphology = collapse_unique(morphology_vme),
    common_id_short = collapse_unique(common_id_short_vme),
    common_id_mid = collapse_unique(common_id_mid_vme),
    common_id_full = collapse_unique(common_id_full_vme),
    .groups = "drop"
  )

write_csv(
  vme_taxonomy,
  file.path(
    directories$vme,
    "vme_type_taxonomy.csv"
  )
)

vme_type_summary <- vme_frame_type_occurrence %>%
  group_by(vme_unit) %>%
  summarise(
    total_vme_annotations = sum(n_vme_annotations),
    frames_present = n_distinct(sample_id_analysis),
    dives_present = n_distinct(dive_id_analysis),
    .groups = "drop"
  ) %>%
  left_join(vme_taxonomy, by = "vme_unit") %>%
  arrange(desc(frames_present), desc(total_vme_annotations), vme_unit)

write_csv(
  vme_type_summary,
  file.path(
    directories$vme,
    "vme_type_occurrence_summary.csv"
  )
)

vme_frame_summary_raw <- vme_data_valid %>%
  group_by(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis
  ) %>%
  summarise(
    n_vme_annotations = n_distinct(point_id_analysis),
    n_vme_types = count_unique(vme_unit),
    vme_types = collapse_unique(vme_unit),
    vme_depth_m = median_or_na(depth_m_vme),
    vme_temperature = median_or_na(temperature_vme),
    .groups = "drop"
  )

write_csv(
  vme_frame_summary_raw,
  file.path(
    directories$vme,
    "vme_frames_observed.csv"
  )
)

vme_dive_summary_raw <- vme_frame_summary_raw %>%
  group_by(dive_id_analysis) %>%
  summarise(
    frames_with_vme = n_distinct(sample_id_analysis),
    total_vme_annotations = sum(n_vme_annotations),
    observed_vme_types = collapse_unique(vme_types),
    .groups = "drop"
  )

write_csv(
  vme_dive_summary_raw,
  file.path(
    directories$vme,
    "vme_occurrence_by_dive_raw.csv"
  )
)

# =============================================================================
# 9. FRAME METADATA, COMPLETENESS, AND COMPOSITION
# =============================================================================

frame_summary <- point_data %>%
  group_by(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis
  ) %>%
  summarise(
    n_points = n(),
    n_biotic = sum(
      point_class == "Biotic"
    ),
    n_abiotic = sum(
      point_class == "Abiotic"
    ),
    n_exclude = sum(
      point_class == "Exclude"
    ),
    n_unsure = sum(
      point_class == "Unsure"
    ),
    n_conflict = sum(
      point_class == "Conflict"
    ),
    n_other = sum(
      point_class == "Other"
    ),
    n_biotic_assigned = sum(
      point_class == "Biotic" &
        !is.na(community_unit) &
        !multiple_community_units
    ),
    n_biotic_unassigned = sum(
      point_class == "Biotic" &
        (
          is.na(community_unit) |
            multiple_community_units
        )
    ),
    depth_range_within_frame_m = range_or_na(
      depth_m
    ),
    temperature_range_within_frame = range_or_na(
      temperature
    ),
    depth_raw_m = median_or_na(
      depth_raw_m
    ),
    depth_m = median_or_na(
      depth_m
    ),
    temperature = median_or_na(
      temperature
    ),
    latitude = median_or_na(
      latitude
    ),
    longitude = median_or_na(
      longitude
    ),
    time = first_non_missing(
      time
    ),
    video_time = first_non_missing(
      video_time
    ),
    .groups = "drop"
  ) %>%
  mutate(
    target_points = TARGET_POINTS_PER_FRAME,
    annotation_completion_pct =
      100 * n_points / target_points,
    below_target = n_points < target_points,
    above_target = n_points > target_points,
    n_valid = n_biotic + n_abiotic,
    pct_biotic_all = if_else(
      n_points > 0,
      100 * n_biotic / n_points,
      NA_real_
    ),
    pct_abiotic_all = if_else(
      n_points > 0,
      100 * n_abiotic / n_points,
      NA_real_
    ),
    pct_exclude_all = if_else(
      n_points > 0,
      100 * n_exclude / n_points,
      NA_real_
    ),
    pct_unsure_all = if_else(
      n_points > 0,
      100 * n_unsure / n_points,
      NA_real_
    ),
    pct_conflict_all = if_else(
      n_points > 0,
      100 * n_conflict / n_points,
      NA_real_
    ),
    pct_biotic_valid = if_else(
      n_valid > 0,
      100 * n_biotic / n_valid,
      NA_real_
    ),
    pct_abiotic_valid = if_else(
      n_valid > 0,
      100 * n_abiotic / n_valid,
      NA_real_
    ),
    pct_biotic_units_resolved = if_else(
      n_biotic > 0,
      100 * n_biotic_assigned / n_biotic,
      NA_real_
    )
  )

# Explicit 100-point QC for the random-point dataset.
frames_not_target_point_count <- frame_summary %>%
  filter(n_points != TARGET_POINTS_PER_FRAME)

write_csv(
  frames_not_target_point_count,
  file.path(
    directories$qc,
    "frames_not_target_point_count.csv"
  )
)

# Join VME occurrence to the complete set of point-data frames. Because the
# point file defines the analysed frame universe, no VME record for a point
# frame is interpreted as no recorded VME observation in that frame.
point_frame_keys <- frame_summary %>%
  select(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis
  )

vme_frames_not_in_point_data <- vme_frame_summary_raw %>%
  anti_join(
    point_frame_keys,
    by = c(
      "dive_id_analysis",
      "frame_id_analysis",
      "sample_id_analysis"
    )
  )

write_csv(
  vme_frames_not_in_point_data,
  file.path(
    directories$vme,
    "vme_frames_not_in_point_data.csv"
  )
)

frame_summary <- frame_summary %>%
  left_join(
    vme_frame_summary_raw %>%
      select(
        dive_id_analysis,
        frame_id_analysis,
        sample_id_analysis,
        n_vme_annotations,
        n_vme_types,
        vme_types
      ),
    by = c(
      "dive_id_analysis",
      "frame_id_analysis",
      "sample_id_analysis"
    )
  ) %>%
  mutate(
    n_vme_annotations = replace_na(n_vme_annotations, 0L),
    n_vme_types = replace_na(n_vme_types, 0L),
    has_vme = n_vme_annotations > 0
  )

write_csv(
  frame_summary %>%
    select(
      dive_id_analysis,
      frame_id_analysis,
      sample_id_analysis,
      has_vme,
      n_vme_annotations,
      n_vme_types,
      vme_types,
      depth_m,
      temperature,
      latitude,
      longitude,
      video_time
    ),
  file.path(
    directories$vme,
    "vme_presence_across_all_point_frames.csv"
  )
)

vme_presence_by_dive <- frame_summary %>%
  group_by(dive_id_analysis) %>%
  summarise(
    n_point_frames = n(),
    frames_with_vme = sum(has_vme),
    frames_without_vme = sum(!has_vme),
    pct_point_frames_with_vme = 100 * frames_with_vme / n_point_frames,
    total_vme_annotations = sum(n_vme_annotations),
    .groups = "drop"
  )

write_csv(
  vme_presence_by_dive,
  file.path(
    directories$vme,
    "vme_presence_by_dive.csv"
  )
)

if (nrow(vme_presence_by_dive) > 0) {
  vme_dive_plot <- ggplot(
    vme_presence_by_dive,
    aes(
      x = dive_id_analysis,
      y = pct_point_frames_with_vme
    )
  ) +
    geom_col() +
    labs(
      title = "Frames with recorded VME observations by dive",
      subtitle = "Denominator = frames represented in point_annotations.csv",
      x = "Dive ID",
      y = "Point-data frames with VME observation (%)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_blank()
    )

  save_plot(
    vme_dive_plot,
    "vme_presence_by_dive",
    directories$vme,
    width = 14,
    height = 8
  )
}

make_vme_presence_plot <- function(data, predictor, output_stub, x_label) {
  plot_data <- data %>%
    filter(!is.na(.data[[predictor]])) %>%
    mutate(vme_present_numeric = as.integer(has_vme))

  if (
    nrow(plot_data) < 10 ||
      n_distinct(plot_data[[predictor]]) < 3 ||
      n_distinct(plot_data$vme_present_numeric) < 2
  ) {
    return(invisible(NULL))
  }

  p <- ggplot(
    plot_data,
    aes(
      x = .data[[predictor]],
      y = vme_present_numeric
    )
  ) +
    geom_jitter(height = 0.04, width = 0, alpha = 0.45) +
    geom_smooth(
      method = "glm",
      method.args = list(family = binomial()),
      se = TRUE
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("No VME recorded", "VME recorded"),
      limits = c(-0.08, 1.08)
    ) +
    labs(
      title = paste0("Recorded VME presence across ", tolower(x_label)),
      subtitle = "Exploratory frame-level logistic smoother; point-data frames define presence/absence universe",
      x = x_label,
      y = "VME observation"
    ) +
    theme_bw(base_size = 11)

  save_plot(
    p,
    output_stub,
    directories$vme,
    width = 9,
    height = 7
  )
}

make_vme_presence_plot(
  frame_summary,
  "depth_m",
  "vme_presence_vs_depth",
  "Depth (m; positive downward)"
)

make_vme_presence_plot(
  frame_summary,
  "temperature",
  "vme_presence_vs_temperature",
  "Temperature"
)

write_csv(
  frame_summary,
  file.path(
    directories$qc,
    "frame_metadata_and_completeness.csv"
  )
)

frame_class_long <- frame_summary %>%
  select(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    starts_with("n_")
  ) %>%
  select(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    n_biotic,
    n_abiotic,
    n_exclude,
    n_unsure,
    n_conflict,
    n_other
  ) %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "point_class",
    values_to = "count"
  ) %>%
  mutate(
    point_class = recode(
      point_class,
      n_biotic = "Biotic",
      n_abiotic = "Abiotic",
      n_exclude = "Exclude",
      n_unsure = "UNSURE",
      n_conflict = "Conflict",
      n_other = "Other"
    )
  )

write_csv(
  frame_class_long,
  file.path(
    directories$composition,
    "frame_point_class_counts.csv"
  )
)

frame_composition_plot <- frame_class_long %>%
  filter(
    point_class %in% c(
      "Biotic",
      "Abiotic",
      "Exclude",
      "UNSURE"
    )
  ) %>%
  group_by(
    dive_id_analysis,
    frame_id_analysis
  ) %>%
  mutate(
    percentage = 100 * count / sum(count)
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = frame_id_analysis,
      y = percentage,
      fill = point_class
    )
  ) +
  geom_col(width = 1) +
  facet_wrap(
    vars(dive_id_analysis),
    scales = "free_x",
    ncol = 4
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Frame composition",
    subtitle = "Biotic, Abiotic, Exclude, and UNSURE point percentages",
    x = "Frame filename",
    y = "Percentage of resolved points",
    fill = "Point class"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

save_plot(
  frame_composition_plot,
  "frame_composition",
  directories$composition,
  width = 22,
  height = max(
    8,
    ceiling(
      n_distinct(frame_summary$dive_id_analysis) / 4
    ) * 3.5
  ),
  dpi = 250
)


# =============================================================================
# 10. FRAME-LEVEL SUBSTRATE AND RELIEF SUMMARIES
# =============================================================================

dominant_value_table <- function(data, value_column, output_name) {
  counts <- data %>%
    filter(
      point_class == "Abiotic",
      !is.na(.data[[value_column]])
    ) %>%
    count(
      dive_id_analysis,
      frame_id_analysis,
      sample_id_analysis,
      value = .data[[value_column]],
      name = "n_points"
    )

  if (nrow(counts) == 0) {
    return(
      tibble(
        sample_id_analysis = character(),
        dominant_value = character(),
        dominant_value_points = integer(),
        dominant_value_pct = numeric()
      )
    )
  }

  totals <- counts %>%
    group_by(sample_id_analysis) %>%
    summarise(
      total_classified_points = sum(n_points),
      .groups = "drop"
    )

  dominant <- counts %>%
    group_by(sample_id_analysis) %>%
    arrange(
      desc(n_points),
      value,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    left_join(
      totals,
      by = "sample_id_analysis"
    ) %>%
    transmute(
      sample_id_analysis,
      dominant_value = value,
      dominant_value_points = n_points,
      dominant_value_pct =
        100 * n_points / total_classified_points
    )

  write_csv(
    counts,
    file.path(
      directories$composition,
      output_name
    )
  )

  dominant
}


dominant_substrate_type <- dominant_value_table(
  point_data,
  "substrate_type",
  "frame_substrate_type_counts.csv"
) %>%
  rename(
    dominant_substrate_type = dominant_value,
    dominant_substrate_type_points =
      dominant_value_points,
    dominant_substrate_type_pct =
      dominant_value_pct
  )

dominant_substrate_detail <- dominant_value_table(
  point_data,
  "substrate_detail",
  "frame_substrate_detail_counts.csv"
) %>%
  rename(
    dominant_substrate_detail = dominant_value,
    dominant_substrate_detail_points =
      dominant_value_points,
    dominant_substrate_detail_pct =
      dominant_value_pct
  )

dominant_substrate_method <- dominant_value_table(
  point_data,
  "substrate_method",
  "frame_substrate_method_counts.csv"
) %>%
  rename(
    dominant_substrate_method = dominant_value,
    dominant_substrate_method_points =
      dominant_value_points,
    dominant_substrate_method_pct =
      dominant_value_pct
  )

dominant_substrate_size <- dominant_value_table(
  point_data,
  "substrate_size",
  "frame_substrate_size_counts.csv"
) %>%
  rename(
    dominant_substrate_size = dominant_value,
    dominant_substrate_size_points =
      dominant_value_points,
    dominant_substrate_size_pct =
      dominant_value_pct
  )

dominant_relief <- dominant_value_table(
  point_data,
  "relief",
  "frame_relief_counts.csv"
) %>%
  rename(
    dominant_relief = dominant_value,
    dominant_relief_points =
      dominant_value_points,
    dominant_relief_pct =
      dominant_value_pct
  )

frame_summary <- frame_summary %>%
  left_join(
    dominant_substrate_type,
    by = "sample_id_analysis"
  ) %>%
  left_join(
    dominant_substrate_detail,
    by = "sample_id_analysis"
  ) %>%
  left_join(
    dominant_substrate_method,
    by = "sample_id_analysis"
  ) %>%
  left_join(
    dominant_substrate_size,
    by = "sample_id_analysis"
  ) %>%
  left_join(
    dominant_relief,
    by = "sample_id_analysis"
  )

write_csv(
  frame_summary,
  file.path(
    directories$composition,
    "frame_ecological_metadata.csv"
  )
)


# =============================================================================
# 11. COMMUNITY MATRIX USING STANDARDISED COMMON IDs
# =============================================================================

community_points <- point_data %>%
  filter(
    point_class == "Biotic",
    !is.na(community_unit),
    !multiple_community_units
  )

if (nrow(community_points) == 0) {
  write_text(
    c(
      "No biological points contained a usable community unit.",
      paste0("Community-unit mode: ", community_unit_mode)
    ),
    file.path(
      directories$logs,
      "pipeline_stopped_no_community_units.txt"
    )
  )

  stop(
    "No usable biological community units were found.",
    call. = FALSE
  )
}

frame_unit_counts <- community_points %>%
  count(
    dive_id_analysis,
    frame_id_analysis,
    sample_id_analysis,
    community_unit,
    name = "count"
  ) %>%
  arrange(
    dive_id_analysis,
    frame_id_analysis,
    community_unit
  )

write_csv(
  frame_unit_counts,
  file.path(
    directories$community,
    "frame_community_unit_counts_long.csv"
  )
)

frame_unit_wide <- frame_unit_counts %>%
  select(
    sample_id_analysis,
    community_unit,
    count
  ) %>%
  pivot_wider(
    names_from = community_unit,
    values_from = count,
    values_fill = 0
  )

write_csv(
  frame_unit_wide,
  file.path(
    directories$community,
    "frame_community_matrix_counts.csv"
  )
)

frame_count_matrix <- as.matrix(
  frame_unit_wide[
    ,
    setdiff(
      names(frame_unit_wide),
      "sample_id_analysis"
    )
  ]
)

rownames(frame_count_matrix) <-
  frame_unit_wide$sample_id_analysis

frame_relative_matrix <- decostand(
  frame_count_matrix,
  method = "total",
  MARGIN = 1
)

write_csv(
  data.frame(
    sample_id_analysis =
      rownames(frame_relative_matrix),
    frame_relative_matrix,
    check.names = FALSE
  ),
  file.path(
    directories$community,
    "frame_community_matrix_relative_abundance.csv"
  )
)


# =============================================================================
# 12. ALPHA DIVERSITY
# =============================================================================

biotic_sample_size <- rowSums(
  frame_count_matrix
)

observed_richness <- specnumber(
  frame_count_matrix
)

shannon <- diversity(
  frame_count_matrix,
  index = "shannon"
)

inverse_simpson <- diversity(
  frame_count_matrix,
  index = "invsimpson"
)

rarefied_richness <- rep(
  NA_real_,
  nrow(frame_count_matrix)
)

eligible_rarefaction <- biotic_sample_size >=
  RAREFACTION_N

if (any(eligible_rarefaction)) {
  rarefied_richness[eligible_rarefaction] <-
    rarefy(
      frame_count_matrix[
        eligible_rarefaction,
        ,
        drop = FALSE
      ],
      sample = RAREFACTION_N
    )
}

singleton_count <- rowSums(
  frame_count_matrix == 1
)

goods_coverage_proxy <- ifelse(
  biotic_sample_size > 0,
  1 - singleton_count / biotic_sample_size,
  NA_real_
)

frame_alpha <- tibble(
  sample_id_analysis =
    rownames(frame_count_matrix),
  assigned_biotic_points =
    biotic_sample_size,
  common_id_richness =
    observed_richness,
  rarefied_common_id_richness =
    rarefied_richness,
  shannon_common_id_diversity =
    shannon,
  inverse_simpson_common_id_diversity =
    inverse_simpson,
  pielou_common_id_evenness = ifelse(
    observed_richness > 1,
    shannon / log(observed_richness),
    NA_real_
  ),
  effective_common_ids_q1 = exp(shannon),
  singleton_common_ids = singleton_count,
  goods_coverage_proxy =
    goods_coverage_proxy,
  eligible_min_biotic_points =
    biotic_sample_size >=
    MIN_BIOTIC_POINTS
) %>%
  left_join(
    frame_summary,
    by = "sample_id_analysis"
  )

write_csv(
  frame_alpha,
  file.path(
    directories$diversity,
    "frame_alpha_diversity.csv"
  )
)

# Retain every frame, including frames with zero assigned biological units.
# Diversity metrics remain NA where they cannot be calculated.
all_frame_metrics <- frame_summary %>%
  left_join(
    frame_alpha %>%
      select(
        sample_id_analysis,
        assigned_biotic_points,
        common_id_richness,
        rarefied_common_id_richness,
        shannon_common_id_diversity,
        inverse_simpson_common_id_diversity,
        pielou_common_id_evenness,
        effective_common_ids_q1,
        singleton_common_ids,
        goods_coverage_proxy,
        eligible_min_biotic_points
      ),
    by = "sample_id_analysis"
  ) %>%
  mutate(
    assigned_biotic_points = replace_na(
      assigned_biotic_points,
      0
    ),
    common_id_richness = replace_na(
      common_id_richness,
      0
    ),
    eligible_min_biotic_points = replace_na(
      eligible_min_biotic_points,
      FALSE
    )
  )

write_csv(
  all_frame_metrics,
  file.path(
    directories$diversity,
    "all_frame_ecological_metrics.csv"
  )
)


# =============================================================================
# 13. DIVE-LEVEL SUMMARY AND HOTSPOT RANKING
# =============================================================================

dive_unit_counts <- community_points %>%
  count(
    dive_id_analysis,
    community_unit,
    name = "count"
  )

dive_unit_wide <- dive_unit_counts %>%
  pivot_wider(
    names_from = community_unit,
    values_from = count,
    values_fill = 0
  )

dive_count_matrix <- as.matrix(
  dive_unit_wide[
    ,
    setdiff(
      names(dive_unit_wide),
      "dive_id_analysis"
    )
  ]
)

rownames(dive_count_matrix) <-
  dive_unit_wide$dive_id_analysis

dive_alpha <- tibble(
  dive_id_analysis =
    rownames(dive_count_matrix),
  assigned_biotic_points =
    rowSums(dive_count_matrix),
  gamma_common_id_richness =
    specnumber(dive_count_matrix),
  shannon_common_id_diversity =
    diversity(
      dive_count_matrix,
      index = "shannon"
    ),
  inverse_simpson_common_id_diversity =
    diversity(
      dive_count_matrix,
      index = "invsimpson"
    )
)

dive_frame_summary <- all_frame_metrics %>%
  group_by(dive_id_analysis) %>%
  summarise(
    n_frames = n(),
    n_frames_with_environment =
      sum(
        !is.na(depth_m) &
          !is.na(temperature)
      ),
    frames_with_vme = sum(has_vme),
    pct_frames_with_vme = 100 * mean(has_vme),
    total_vme_annotations = sum(n_vme_annotations),
    mean_vme_types_per_frame = mean(n_vme_types),
    mean_biotic_cover_valid =
      mean(
        pct_biotic_valid,
        na.rm = TRUE
      ),
    median_biotic_cover_valid =
      median(
        pct_biotic_valid,
        na.rm = TRUE
      ),
    mean_frame_common_id_richness =
      mean(
        common_id_richness,
        na.rm = TRUE
      ),
    mean_frame_rarefied_common_id_richness =
      mean(
        rarefied_common_id_richness,
        na.rm = TRUE
      ),
    mean_frame_shannon_common_id =
      mean(
        shannon_common_id_diversity,
        na.rm = TRUE
      ),
    mean_exclude_pct =
      mean(
        pct_exclude_all,
        na.rm = TRUE
      ),
    min_depth_m = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      min(depth_m, na.rm = TRUE)
    ),
    max_depth_m = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      max(depth_m, na.rm = TRUE)
    ),
    min_temperature = ifelse(
      all(is.na(temperature)),
      NA_real_,
      min(temperature, na.rm = TRUE)
    ),
    max_temperature = ifelse(
      all(is.na(temperature)),
      NA_real_,
      max(temperature, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  left_join(
    dive_alpha,
    by = "dive_id_analysis"
  )

write_csv(
  dive_frame_summary,
  file.path(
    directories$diversity,
    "dive_ecological_summary.csv"
  )
)

hotspot_frames <- all_frame_metrics %>%
  filter(
    eligible_min_biotic_points
  ) %>%
  mutate(
    richness_z = scale_safe(
      rarefied_common_id_richness
    ),
    shannon_z = scale_safe(
      shannon_common_id_diversity
    ),
    biotic_cover_z = scale_safe(
      pct_biotic_valid
    ),
    hotspot_score = rowMeans(
      cbind(
        richness_z,
        shannon_z,
        biotic_cover_z
      ),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(hotspot_score))

write_csv(
  hotspot_frames,
  file.path(
    directories$diversity,
    "candidate_frame_hotspot_ranking.csv"
  )
)


# =============================================================================
# 14. BASIC DIVERSITY PLOTS
# =============================================================================

richness_by_dive_plot <- frame_alpha %>%
  filter(
    eligible_min_biotic_points
  ) %>%
  ggplot(
    aes(
      x = dive_id_analysis,
      y = rarefied_common_id_richness
    )
  ) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.55,
    size = 1.5
  ) +
  labs(
    title = "Standardised common-ID richness by dive",
    subtitle = paste0(
      "Rarefied to ",
      RAREFACTION_N,
      " assigned biological points per frame"
    ),
    x = "Dive ID",
    y = "Rarefied group richness"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.major.x = element_blank()
  )

save_plot(
  richness_by_dive_plot,
  "rarefied_common_id_richness_by_dive",
  directories$diversity,
  width = 15,
  height = 8
)

shannon_by_dive_plot <- frame_alpha %>%
  filter(
    eligible_min_biotic_points
  ) %>%
  ggplot(
    aes(
      x = dive_id_analysis,
      y = shannon_common_id_diversity
    )
  ) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.55,
    size = 1.5
  ) +
  labs(
    title = "Common-ID group diversity by dive",
    x = "Dive ID",
    y = "Shannon group diversity"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.major.x = element_blank()
  )

save_plot(
  shannon_by_dive_plot,
  "shannon_common_id_diversity_by_dive",
  directories$diversity,
  width = 15,
  height = 8
)


# =============================================================================
# 15. ENVIRONMENTAL COVERAGE AND CONFOUNDING CHECKS
# =============================================================================

environment_coverage <- all_frame_metrics %>%
  group_by(dive_id_analysis) %>%
  summarise(
    n_frames = n(),
    n_depth = sum(!is.na(depth_m)),
    n_temperature =
      sum(!is.na(temperature)),
    depth_min = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      min(depth_m, na.rm = TRUE)
    ),
    depth_max = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      max(depth_m, na.rm = TRUE)
    ),
    temperature_min = ifelse(
      all(is.na(temperature)),
      NA_real_,
      min(temperature, na.rm = TRUE)
    ),
    temperature_max = ifelse(
      all(is.na(temperature)),
      NA_real_,
      max(temperature, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write_csv(
  environment_coverage,
  file.path(
    directories$environment,
    "environmental_data_coverage_by_dive.csv"
  )
)

complete_environment <- all_frame_metrics %>%
  filter(
    !is.na(depth_m),
    !is.na(temperature)
  )

depth_temperature_cor <- safe_spearman(
  complete_environment$depth_m,
  complete_environment$temperature
)

write_csv(
  depth_temperature_cor,
  file.path(
    directories$environment,
    "depth_temperature_spearman_correlation.csv"
  )
)


# Q1: Is observation quality confounded with the environment?

make_environment_plot(
  all_frame_metrics,
  "depth_m",
  "pct_exclude_all",
  "Does excluded image area change with depth?",
  "Exclude points (%)",
  "exclude_percentage_vs_depth",
  directories$environment,
  percent_axis = TRUE
)

make_environment_plot(
  all_frame_metrics,
  "temperature",
  "pct_exclude_all",
  "Does excluded image area change with temperature?",
  "Exclude points (%)",
  "exclude_percentage_vs_temperature",
  directories$environment,
  percent_axis = TRUE
)


# Q2: Does biotic cover change with the environment?

make_environment_plot(
  all_frame_metrics,
  "depth_m",
  "pct_biotic_valid",
  "Does biotic cover change with depth?",
  "Biotic percentage of valid points",
  "biotic_cover_vs_depth",
  directories$environment,
  percent_axis = TRUE
)

make_environment_plot(
  all_frame_metrics,
  "temperature",
  "pct_biotic_valid",
  "Does biotic cover change with temperature?",
  "Biotic percentage of valid points",
  "biotic_cover_vs_temperature",
  directories$environment,
  percent_axis = TRUE
)


# Q3: Does group richness/diversity change with the environment?

eligible_alpha <- all_frame_metrics %>%
  filter(
    eligible_min_biotic_points
  )

make_environment_plot(
  eligible_alpha,
  "depth_m",
  "rarefied_common_id_richness",
  "Does standardised group richness change with depth?",
  paste0(
    "Rarefied group richness (n = ",
    RAREFACTION_N,
    ")"
  ),
  "rarefied_richness_vs_depth",
  directories$environment
)

make_environment_plot(
  eligible_alpha,
  "temperature",
  "rarefied_common_id_richness",
  "Does standardised group richness change with temperature?",
  paste0(
    "Rarefied group richness (n = ",
    RAREFACTION_N,
    ")"
  ),
  "rarefied_richness_vs_temperature",
  directories$environment
)

make_environment_plot(
  eligible_alpha,
  "depth_m",
  "shannon_common_id_diversity",
  "Does common-ID group diversity change with depth?",
  "Shannon group diversity",
  "shannon_diversity_vs_depth",
  directories$environment
)

make_environment_plot(
  eligible_alpha,
  "temperature",
  "shannon_common_id_diversity",
  "Does common-ID group diversity change with temperature?",
  "Shannon group diversity",
  "shannon_diversity_vs_temperature",
  directories$environment
)


# Q2b: Does biotic cover vary among dominant substrate types?

substrate_cover_data <- all_frame_metrics %>%
  filter(
    !is.na(dominant_substrate_type),
    !is.na(pct_biotic_valid)
  )

if (
  nrow(substrate_cover_data) >= 6 &&
    n_distinct(
      substrate_cover_data$dominant_substrate_type
    ) >= 2
) {
  substrate_cover_plot <- ggplot(
    substrate_cover_data,
    aes(
      x = dominant_substrate_type,
      y = pct_biotic_valid
    )
  ) +
    geom_boxplot(
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.65
    ) +
    labs(
      title = "Biotic cover among dominant substrate types",
      subtitle = "Descriptive association; substrate and cover derive from the same point set",
      x = "Dominant substrate type",
      y = "Biotic percentage of valid points"
    ) +
    theme_bw(base_size = 11)

  save_plot(
    substrate_cover_plot,
    "biotic_cover_by_dominant_substrate",
    directories$environment,
    width = 8,
    height = 7
  )
}


# =============================================================================
# 16. GROUP OCCUPANCY AND ENVIRONMENTAL ASSOCIATIONS
# =============================================================================

group_frame_data <- frame_unit_counts %>%
  left_join(
    frame_summary %>%
      select(
        sample_id_analysis,
        dive_id_analysis,
        frame_id_analysis,
        n_biotic_assigned,
        depth_m,
        temperature,
        dominant_substrate_type,
        dominant_substrate_detail
      ),
    by = c(
      "sample_id_analysis",
      "dive_id_analysis",
      "frame_id_analysis"
    )
  ) %>%
  mutate(
    relative_abundance = if_else(
      n_biotic_assigned > 0,
      count / n_biotic_assigned,
      NA_real_
    )
  )

group_taxonomy <- community_points %>%
  group_by(community_unit) %>%
  summarise(
    cpc_codes = collapse_unique(cpc_codes),
    kingdom = collapse_unique(kingdom),
    phylum = collapse_unique(phylum),
    class = collapse_unique(taxon_class),
    order = collapse_unique(order),
    family = collapse_unique(family),
    taxonomic_resolution =
      collapse_unique(taxonomic_resolution),
    morphology = collapse_unique(morphology),
    common_id_short = collapse_unique(common_id_short),
    common_id_mid = collapse_unique(common_id_mid),
    common_id_full = collapse_unique(common_id_full),
    .groups = "drop"
  )

write_csv(
  group_taxonomy,
  file.path(
    directories$groups,
    "community_unit_taxonomy.csv"
  )
)

group_summary <- group_frame_data %>%
  group_by(community_unit) %>%
  summarise(
    total_points = sum(count),
    frames_present =
      n_distinct(sample_id_analysis),
    dives_present =
      n_distinct(dive_id_analysis),
    frame_occupancy_pct =
      100 * frames_present /
      n_distinct(frame_summary$sample_id_analysis),
    weighted_mean_depth_m =
      weighted_mean_or_na(
        depth_m,
        count
      ),
    min_depth_m = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      min(depth_m, na.rm = TRUE)
    ),
    max_depth_m = ifelse(
      all(is.na(depth_m)),
      NA_real_,
      max(depth_m, na.rm = TRUE)
    ),
    weighted_mean_temperature =
      weighted_mean_or_na(
        temperature,
        count
      ),
    min_temperature = ifelse(
      all(is.na(temperature)),
      NA_real_,
      min(temperature, na.rm = TRUE)
    ),
    max_temperature = ifelse(
      all(is.na(temperature)),
      NA_real_,
      max(temperature, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  left_join(
    group_taxonomy,
    by = "community_unit"
  ) %>%
  arrange(desc(total_points))

write_csv(
  group_summary,
  file.path(
    directories$groups,
    "community_group_occupancy_and_environment.csv"
  )
)

top_groups <- head(
  group_summary$community_unit,
  TOP_N_GROUPS
)

environment_frames <- frame_summary %>%
  filter(
    n_biotic_assigned > 0,
    !is.na(depth_m) |
      !is.na(temperature)
  ) %>%
  select(
    sample_id_analysis,
    dive_id_analysis,
    frame_id_analysis,
    n_biotic_assigned,
    depth_m,
    temperature
  )

top_group_grid <- expand_grid(
  sample_id_analysis =
    environment_frames$sample_id_analysis,
  community_unit = top_groups
) %>%
  left_join(
    group_frame_data %>%
      filter(
        community_unit %in% top_groups
      ) %>%
      select(
        sample_id_analysis,
        community_unit,
        count
      ),
    by = c(
      "sample_id_analysis",
      "community_unit"
    )
  ) %>%
  mutate(
    count = replace_na(
      count,
      0L
    )
  ) %>%
  left_join(
    environment_frames,
    by = "sample_id_analysis"
  ) %>%
  mutate(
    relative_abundance = if_else(
      n_biotic_assigned > 0,
      count / n_biotic_assigned,
      NA_real_
    )
  )

write_csv(
  top_group_grid,
  file.path(
    directories$groups,
    "top_group_frame_relative_abundance.csv"
  )
)


# Q6: Which groups show depth or temperature associations?

group_correlation_rows <- list()

for (group_name in top_groups) {
  group_data <- top_group_grid %>%
    filter(
      community_unit == group_name
    )

  depth_result <- safe_spearman(
    group_data$relative_abundance,
    group_data$depth_m
  )

  temperature_result <- safe_spearman(
    group_data$relative_abundance,
    group_data$temperature
  )

  group_correlation_rows[[
    paste0(group_name, "_depth")
  ]] <- data.frame(
    community_unit = group_name,
    predictor = "depth_m",
    rho = depth_result$rho,
    p_value = depth_result$p_value,
    n = depth_result$n
  )

  group_correlation_rows[[
    paste0(group_name, "_temperature")
  ]] <- data.frame(
    community_unit = group_name,
    predictor = "temperature",
    rho = temperature_result$rho,
    p_value = temperature_result$p_value,
    n = temperature_result$n
  )
}

group_correlations <- bind_rows(
  group_correlation_rows
) %>%
  group_by(predictor) %>%
  mutate(
    p_adjusted_bh = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  arrange(
    predictor,
    p_adjusted_bh,
    desc(abs(rho))
  )

write_csv(
  group_correlations,
  file.path(
    directories$groups,
    "exploratory_group_environment_correlations.csv"
  )
)

if (
  nrow(top_group_grid) > 0 &&
    any(!is.na(top_group_grid$depth_m))
) {
  group_depth_plot <- top_group_grid %>%
    filter(
      !is.na(depth_m)
    ) %>%
    ggplot(
      aes(
        x = depth_m,
        y = relative_abundance,
        colour = dive_id_analysis
      )
    ) +
    geom_point(
      alpha = 0.65,
      size = 1.3
    ) +
    facet_wrap(
      vars(community_unit),
      scales = "free_y",
      ncol = 4
    ) +
    labs(
      title = "Common biological groups across depth",
      subtitle = "Relative abundance among assigned biological points",
      x = "Depth (m; positive downward)",
      y = "Relative abundance",
      colour = "Dive ID"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )

  save_plot(
    group_depth_plot,
    "top_group_relative_abundance_vs_depth",
    directories$groups,
    width = 16,
    height = 12
  )
}

if (
  nrow(top_group_grid) > 0 &&
    any(!is.na(top_group_grid$temperature))
) {
  group_temperature_plot <- top_group_grid %>%
    filter(
      !is.na(temperature)
    ) %>%
    ggplot(
      aes(
        x = temperature,
        y = relative_abundance,
        colour = dive_id_analysis
      )
    ) +
    geom_point(
      alpha = 0.65,
      size = 1.3
    ) +
    facet_wrap(
      vars(community_unit),
      scales = "free_y",
      ncol = 4
    ) +
    labs(
      title = "Common biological groups across temperature",
      subtitle = "Relative abundance among assigned biological points",
      x = "Temperature",
      y = "Relative abundance",
      colour = "Dive ID"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )

  save_plot(
    group_temperature_plot,
    "top_group_relative_abundance_vs_temperature",
    directories$groups,
    width = 16,
    height = 12
  )
}


# Frame-by-group heatmap.

frame_order_table <- frame_summary %>%
  arrange(
    dive_id_analysis,
    coalesce(
      depth_m,
      Inf
    ),
    frame_id_analysis
  ) %>%
  mutate(
    frame_order = row_number()
  )

heatmap_data <- expand_grid(
  sample_id_analysis =
    frame_order_table$sample_id_analysis,
  community_unit = rev(top_groups)
) %>%
  left_join(
    group_frame_data %>%
      filter(
        community_unit %in% top_groups
      ) %>%
      select(
        sample_id_analysis,
        community_unit,
        relative_abundance
      ),
    by = c(
      "sample_id_analysis",
      "community_unit"
    )
  ) %>%
  mutate(
    relative_abundance = replace_na(
      relative_abundance,
      0
    )
  ) %>%
  left_join(
    frame_order_table %>%
      select(
        sample_id_analysis,
        dive_id_analysis,
        frame_order
      ),
    by = "sample_id_analysis"
  )

if (nrow(heatmap_data) > 0) {
  group_heatmap <- ggplot(
    heatmap_data,
    aes(
      x = frame_order,
      y = factor(
        community_unit,
        levels = rev(top_groups)
      ),
      fill = relative_abundance
    )
  ) +
    geom_tile() +
    facet_grid(
      cols = vars(dive_id_analysis),
      scales = "free_x",
      space = "free_x"
    ) +
    labs(
      title = "Relative abundance of common biological groups",
      subtitle = "Frames are ordered by dive and depth where available",
      x = "Frame order",
      y = "Community group",
      fill = "Relative\nabundance"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(
        angle = 90,
        hjust = 0
      )
    )

  save_plot(
    group_heatmap,
    "top_group_relative_abundance_heatmap",
    directories$groups,
    width = 20,
    height = 8
  )
}


# =============================================================================
# 17. BRAY-CURTIS NMDS AND COMMUNITY TURNOVER
# =============================================================================

ordination_keep <- rowSums(
  frame_count_matrix
) >= MIN_BIOTIC_POINTS

ordination_matrix <- frame_relative_matrix[
  ordination_keep,
  ,
  drop = FALSE
]

ordination_matrix <- ordination_matrix[
  ,
  colSums(ordination_matrix) > 0,
  drop = FALSE
]

nmds <- NULL
nmds_scores <- NULL
bray_distance <- NULL
nmds_stress <- NA_real_
envfit_result <- NULL
permanova_results <- list()

if (
  nrow(ordination_matrix) >= 4 &&
    ncol(ordination_matrix) >= 2
) {
  bray_distance <- vegdist(
    ordination_matrix,
    method = "bray"
  )

  nmds <- tryCatch(
    metaMDS(
      ordination_matrix,
      distance = "bray",
      k = 2,
      trymax = 200,
      autotransform = FALSE,
      trace = FALSE
    ),
    error = function(error) {
      warning(
        paste0(
          "NMDS failed: ",
          error$message
        ),
        call. = FALSE
      )

      NULL
    }
  )

  if (!is.null(nmds)) {
    nmds_stress <- nmds$stress

    nmds_scores <- as.data.frame(
      scores(
        nmds,
        display = "sites"
      )
    )

    nmds_scores$sample_id_analysis <-
      rownames(nmds_scores)

    nmds_scores <- nmds_scores %>%
      left_join(
        frame_summary,
        by = "sample_id_analysis"
      )

    write_csv(
      nmds_scores,
      file.path(
        directories$ordination,
        "frame_nmds_scores.csv"
      )
    )

    write_text(
      paste0(
        "NMDS stress: ",
        round(nmds_stress, 4)
      ),
      file.path(
        directories$ordination,
        "frame_nmds_stress.txt"
      )
    )

    nmds_dive_plot <- ggplot(
      nmds_scores,
      aes(
        x = NMDS1,
        y = NMDS2,
        colour = dive_id_analysis
      )
    ) +
      geom_point(
        size = 2.4,
        alpha = 0.8
      ) +
      coord_equal() +
      labs(
        title = "Frame-level assemblage composition",
        subtitle = paste0(
          "Bray-Curtis NMDS; stress = ",
          round(nmds_stress, 3)
        ),
        colour = "Dive ID"
      ) +
      theme_bw(base_size = 11)

    save_plot(
      nmds_dive_plot,
      "frame_nmds_by_dive",
      directories$ordination,
      width = 10,
      height = 8
    )

    if (any(!is.na(nmds_scores$depth_m))) {
      nmds_depth_plot <- ggplot(
        nmds_scores %>%
          filter(
            !is.na(depth_m)
          ),
        aes(
          x = NMDS1,
          y = NMDS2,
          colour = depth_m
        )
      ) +
        geom_point(
          size = 2.6,
          alpha = 0.85
        ) +
        coord_equal() +
        labs(
          title = "Assemblage composition across depth",
          subtitle = paste0(
            "Bray-Curtis NMDS; stress = ",
            round(nmds_stress, 3)
          ),
          colour = "Depth (m)"
        ) +
        theme_bw(base_size = 11)

      save_plot(
        nmds_depth_plot,
        "frame_nmds_by_depth",
        directories$ordination,
        width = 10,
        height = 8
      )
    }

    if (any(!is.na(nmds_scores$temperature))) {
      nmds_temperature_plot <- ggplot(
        nmds_scores %>%
          filter(
            !is.na(temperature)
          ),
        aes(
          x = NMDS1,
          y = NMDS2,
          colour = temperature
        )
      ) +
        geom_point(
          size = 2.6,
          alpha = 0.85
        ) +
        coord_equal() +
        labs(
          title = "Assemblage composition across temperature",
          subtitle = paste0(
            "Bray-Curtis NMDS; stress = ",
            round(nmds_stress, 3)
          ),
          colour = "Temperature"
        ) +
        theme_bw(base_size = 11)

      save_plot(
        nmds_temperature_plot,
        "frame_nmds_by_temperature",
        directories$ordination,
        width = 10,
        height = 8
      )
    }

    substrate_nmds_data <- nmds_scores %>%
      filter(
        !is.na(dominant_substrate_type)
      )

    if (
      nrow(substrate_nmds_data) >= 6 &&
        n_distinct(
          substrate_nmds_data$dominant_substrate_type
        ) >= 2
    ) {
      nmds_substrate_plot <- ggplot(
        substrate_nmds_data,
        aes(
          x = NMDS1,
          y = NMDS2,
          colour = dominant_substrate_type
        )
      ) +
        geom_point(
          size = 2.6,
          alpha = 0.85
        ) +
        coord_equal() +
        labs(
          title = "Assemblage composition by dominant substrate",
          subtitle = paste0(
            "Bray-Curtis NMDS; stress = ",
            round(nmds_stress, 3)
          ),
          colour = "Dominant substrate"
        ) +
        theme_bw(base_size = 11)

      save_plot(
        nmds_substrate_plot,
        "frame_nmds_by_substrate",
        directories$ordination,
        width = 10,
        height = 8
      )
    }

    # Environmental vector fitting.
    envfit_metadata <- nmds_scores %>%
      select(
        sample_id_analysis,
        depth_m,
        temperature
      )

    envfit_keep <- complete.cases(
      envfit_metadata[
        ,
        c(
          "depth_m",
          "temperature"
        )
      ]
    )

    if (sum(envfit_keep) >= 10) {
      envfit_matrix <- ordination_matrix[
        envfit_metadata$sample_id_analysis[
          envfit_keep
        ],
        ,
        drop = FALSE
      ]

      envfit_ordination <- metaMDS(
        envfit_matrix,
        distance = "bray",
        k = 2,
        trymax = 100,
        autotransform = FALSE,
        trace = FALSE
      )

      envfit_result <- envfit(
        envfit_ordination,
        envfit_metadata[
          envfit_keep,
          c(
            "depth_m",
            "temperature"
          )
        ],
        permutations = N_PERMUTATIONS
      )

      capture.output(
        envfit_result,
        file = file.path(
          directories$ordination,
          "nmds_environmental_vector_fit.txt"
        )
      )
    }
  }
}


# =============================================================================
# 18. EXPLORATORY GAM MODELS
# =============================================================================
#
# These models are hypothesis-generating. Adjacent frames may be temporally
# autocorrelated and current environmental metadata may cover only a subset of
# dives. Separate depth and temperature models are fitted to avoid forcing two
# potentially collinear predictors into one model.
# =============================================================================

model_registry <- list()


fit_exploratory_gam <- function(
    data,
    response_type,
    predictor,
    model_name
) {
  model_data <- data %>%
    filter(
      !is.na(.data[[predictor]])
    ) %>%
    mutate(
      dive_factor = factor(
        dive_id_analysis
      )
    )

  if (response_type == "biotic_cover") {
    model_data <- model_data %>%
      filter(
        n_biotic + n_abiotic > 0
      )
  } else if (response_type == "rarefied_richness") {
    model_data <- model_data %>%
      filter(
        !is.na(rarefied_common_id_richness)
      )
  } else if (response_type == "shannon") {
    model_data <- model_data %>%
      filter(
        !is.na(shannon_common_id_diversity)
      )
  } else if (response_type == "exclude") {
    model_data <- model_data %>%
      filter(
        n_points > 0
      )
  }

  if (
    nrow(model_data) < 15 ||
      n_distinct(model_data[[predictor]]) < 6
  ) {
    return(NULL)
  }

  k_value <- min(
    5,
    max(
      3,
      n_distinct(
        model_data[[predictor]]
      ) - 1
    )
  )

  dive_term <- if (
    n_distinct(model_data$dive_factor) >= 3
  ) {
    " + s(dive_factor, bs = 're')"
  } else if (
    n_distinct(model_data$dive_factor) == 2
  ) {
    " + dive_factor"
  } else {
    ""
  }

  if (response_type == "biotic_cover") {
    formula_text <- paste0(
      "cbind(n_biotic, n_abiotic) ~ ",
      "s(",
      predictor,
      ", bs = 'cs', k = ",
      k_value,
      ")",
      dive_term
    )

    model_family <- quasibinomial()
  } else if (response_type == "exclude") {
    formula_text <- paste0(
      "cbind(n_exclude, pmax(n_points - n_exclude, 0)) ~ ",
      "s(",
      predictor,
      ", bs = 'cs', k = ",
      k_value,
      ")",
      dive_term
    )

    model_family <- quasibinomial()
  } else if (
    response_type == "rarefied_richness"
  ) {
    formula_text <- paste0(
      "rarefied_common_id_richness ~ ",
      "s(",
      predictor,
      ", bs = 'cs', k = ",
      k_value,
      ")",
      dive_term
    )

    model_family <- gaussian()
  } else {
    formula_text <- paste0(
      "shannon_common_id_diversity ~ ",
      "s(",
      predictor,
      ", bs = 'cs', k = ",
      k_value,
      ")",
      dive_term
    )

    model_family <- gaussian()
  }

  fitted_model <- tryCatch(
    gam(
      formula = as.formula(
        formula_text
      ),
      data = model_data,
      family = model_family,
      method = "REML"
    ),
    error = function(error) {
      warning(
        paste0(
          "Model ",
          model_name,
          " failed: ",
          error$message
        ),
        call. = FALSE
      )

      NULL
    }
  )

  if (is.null(fitted_model)) {
    return(NULL)
  }

  model_summary <- summary(
    fitted_model
  )

  predictor_p <- if (
    !is.null(model_summary$s.table) &&
      nrow(model_summary$s.table) >= 1
  ) {
    model_summary$s.table[
      1,
      "p-value"
    ]
  } else {
    NA_real_
  }

  capture.output(
    model_summary,
    file = file.path(
      directories$models,
      paste0(
        model_name,
        "_summary.txt"
      )
    )
  )

  saveRDS(
    fitted_model,
    file.path(
      directories$models,
      paste0(
        model_name,
        ".rds"
      )
    )
  )

  data.frame(
    model = model_name,
    response = response_type,
    predictor = predictor,
    n_frames = nrow(model_data),
    n_dives =
      n_distinct(model_data$dive_factor),
    smooth_p_value = predictor_p,
    formula = formula_text,
    stringsAsFactors = FALSE
  )
}


if (RUN_EXPLORATORY_GAMS) {
  gam_specs <- list(
    c(
      "exclude",
      "depth_m",
      "exclude_vs_depth_gam"
    ),
    c(
      "exclude",
      "temperature",
      "exclude_vs_temperature_gam"
    ),
    c(
      "biotic_cover",
      "depth_m",
      "biotic_cover_vs_depth_gam"
    ),
    c(
      "biotic_cover",
      "temperature",
      "biotic_cover_vs_temperature_gam"
    ),
    c(
      "rarefied_richness",
      "depth_m",
      "rarefied_richness_vs_depth_gam"
    ),
    c(
      "rarefied_richness",
      "temperature",
      "rarefied_richness_vs_temperature_gam"
    ),
    c(
      "shannon",
      "depth_m",
      "shannon_vs_depth_gam"
    ),
    c(
      "shannon",
      "temperature",
      "shannon_vs_temperature_gam"
    )
  )

  for (specification in gam_specs) {
    fitted_record <- fit_exploratory_gam(
      all_frame_metrics,
      response_type = specification[[1]],
      predictor = specification[[2]],
      model_name = specification[[3]]
    )

    if (!is.null(fitted_record)) {
      model_registry[[
        specification[[3]]
      ]] <- fitted_record
    }
  }
}

model_registry_table <- bind_rows(
  model_registry
)

write_csv(
  model_registry_table,
  file.path(
    directories$models,
    "exploratory_gam_registry.csv"
  )
)


# =============================================================================
# 19. EXPLORATORY PERMANOVA
# =============================================================================

if (
  RUN_PERMANOVA &&
    !is.null(bray_distance) &&
    nrow(ordination_matrix) >= 8
) {
  permanova_metadata <- frame_summary %>%
    filter(
      sample_id_analysis %in%
        rownames(ordination_matrix)
    ) %>%
    arrange(
      match(
        sample_id_analysis,
        rownames(ordination_matrix)
      )
    )

  run_permanova_predictor <- function(
      predictor,
      output_name
  ) {
    keep <- !is.na(
      permanova_metadata[[predictor]]
    )

    if (
      sum(keep) < 8 ||
        n_distinct(
          permanova_metadata[[predictor]][keep]
        ) < 2
    ) {
      return(NULL)
    }

    matrix_subset <- ordination_matrix[
      keep,
      ,
      drop = FALSE
    ]

    metadata_subset <-
      permanova_metadata[
        keep,
        ,
        drop = FALSE
      ]

    distance_subset <- vegdist(
      matrix_subset,
      method = "bray"
    )

    formula_object <- as.formula(
      paste0(
        "distance_subset ~ ",
        predictor
      )
    )

    if (
      n_distinct(
        metadata_subset$dive_id_analysis
      ) >= 2
    ) {
      result <- adonis2(
        formula_object,
        data = metadata_subset,
        permutations = N_PERMUTATIONS,
        strata =
          metadata_subset$dive_id_analysis
      )
    } else {
      result <- adonis2(
        formula_object,
        data = metadata_subset,
        permutations = N_PERMUTATIONS
      )
    }

    capture.output(
      result,
      file = file.path(
        directories$models,
        paste0(
          output_name,
          ".txt"
        )
      )
    )

    result
  }

  permanova_results$depth <-
    run_permanova_predictor(
      "depth_m",
      "permanova_community_vs_depth"
    )

  permanova_results$temperature <-
    run_permanova_predictor(
      "temperature",
      "permanova_community_vs_temperature"
    )

  permanova_results$substrate <-
    run_permanova_predictor(
      "dominant_substrate_type",
      "permanova_community_vs_substrate"
    )

  substrate_keep <- !is.na(
    permanova_metadata$dominant_substrate_type
  )

  substrate_groups <- table(
    permanova_metadata$dominant_substrate_type[
      substrate_keep
    ]
  )

  if (
    sum(substrate_keep) >= 8 &&
      sum(substrate_groups >= 3) >= 2
  ) {
    substrate_matrix <- ordination_matrix[
      substrate_keep,
      ,
      drop = FALSE
    ]

    substrate_distance <- vegdist(
      substrate_matrix,
      method = "bray"
    )

    dispersion <- betadisper(
      substrate_distance,
      group = factor(
        permanova_metadata$dominant_substrate_type[
          substrate_keep
        ]
      )
    )

    dispersion_test <- permutest(
      dispersion,
      permutations = N_PERMUTATIONS
    )

    capture.output(
      dispersion_test,
      file = file.path(
        directories$models,
        "multivariate_dispersion_by_substrate.txt"
      )
    )
  }
}


# =============================================================================
# 20. AUTOMATED ECOLOGICAL QUESTION SUMMARY
# =============================================================================

top_hotspots <- head(
  hotspot_frames,
  10
)

strong_group_associations <- group_correlations %>%
  filter(
    !is.na(rho)
  ) %>%
  arrange(
    p_adjusted_bh,
    desc(abs(rho))
  ) %>%
  slice_head(n = 10)

summary_lines <- c(
  "BIIGLE ECOLOGICAL ANALYSIS SUMMARY",
  "==================================",
  "",
  paste0(
    "Point-data input: ",
    normalizePath(point_input_csv)
  ),
  paste0(
    "VME input: ",
    normalizePath(vme_input_csv)
  ),
  paste0(
    "Community-unit mode: ",
    community_unit_mode
  ),
  paste0(
    "Rows with common_id_full: ",
    sum(!is.na(common_id_full_raw))
  ),
  paste0(
    "Rows with depth_m: ",
    sum(!is.na(safe_numeric(column_or_na(raw_data, DEPTH_COL, numeric = TRUE))))
  ),
  paste0(
    "Rows with temperature: ",
    sum(!is.na(safe_numeric(column_or_na(raw_data, TEMPERATURE_COL, numeric = TRUE))))
  ),
  paste0(
    "Point-data rows read: ",
    format(nrow(point_raw_data), big.mark = ",")
  ),
  paste0(
    "VME rows read: ",
    format(nrow(vme_raw_data), big.mark = ",")
  ),
  paste0(
    "Unique random annotation points: ",
    format(
      nrow(point_data),
      big.mark = ","
    )
  ),
  paste0(
    "Frames: ",
    n_distinct(
      point_data$sample_id_analysis
    )
  ),
  paste0(
    "Dives: ",
    n_distinct(
      point_data$dive_id_analysis
    )
  ),
  paste0(
    "Point frames not exactly ",
    TARGET_POINTS_PER_FRAME,
    " unique points: ",
    nrow(frames_not_target_point_count)
  ),
  paste0(
    "Frames with depth: ",
    sum(
      !is.na(frame_summary$depth_m)
    )
  ),
  paste0(
    "Frames with temperature: ",
    sum(
      !is.na(frame_summary$temperature)
    )
  ),
  paste0(
    "Point-data frames with recorded VME: ",
    sum(frame_summary$has_vme)
  ),
  paste0(
    "Unique VME annotations: ",
    nrow(vme_annotations_unique)
  ),
  paste0(
    "VME types: ",
    n_distinct(vme_frame_type_occurrence$vme_unit)
  ),
  paste0(
    "VME frames not represented in point data: ",
    nrow(vme_frames_not_in_point_data)
  ),
  paste0(
    "Multiply labelled exported points: ",
    nrow(duplicate_export_points)
  ),
  paste0(
    "Conflicting top-level points excluded: ",
    nrow(conflicting_points)
  ),
  paste0(
    "Biotic points with multiple community units excluded from community matrices: ",
    nrow(multiple_unit_points)
  ),
  "",
  "QUESTIONS AND OUTPUTS",
  "---------------------",
  "",
  "Q1. Is image quality confounded with depth or temperature?",
  "    See 05_environment/exclude_percentage_vs_depth.*",
  "    See 05_environment/exclude_percentage_vs_temperature.*",
  "    See 09_models/exclude_vs_*_gam_summary.txt",
  "",
  "Q2. Does biotic cover vary with depth, temperature, or substrate?",
  "    See 05_environment/biotic_cover_vs_depth.*",
  "    See 05_environment/biotic_cover_vs_temperature.*",
  "    See 05_environment/biotic_cover_by_dominant_substrate.*",
  "    See 09_models/biotic_cover_vs_*_gam_summary.txt",
  "",
  "Q3. Does standardised group richness or Shannon diversity vary environmentally?",
  "    See 05_environment/rarefied_richness_vs_*.*",
  "    See 05_environment/shannon_diversity_vs_*.*",
  "    See 09_models/rarefied_richness_vs_*_gam_summary.txt",
  "    See 09_models/shannon_vs_*_gam_summary.txt",
  "",
  "Q4. Does assemblage composition turn over among frames and dives?",
  "    See 06_ordination/frame_nmds_by_dive.*",
  "    See 06_ordination/frame_nmds_by_depth.*",
  "    See 06_ordination/frame_nmds_by_temperature.*",
  "",
  "Q5. Are depth, temperature, or substrate associated with assemblage composition?",
  "    See 06_ordination/nmds_environmental_vector_fit.txt",
  "    See 09_models/permanova_community_vs_*.txt",
  "",
  "Q6. Which biological groups show depth or temperature associations?",
  "    See 07_group_responses/exploratory_group_environment_correlations.csv",
  "    See 07_group_responses/community_group_occupancy_and_environment.csv",
  "",
  "Q7. Which frames are candidate biodiversity hotspots?",
  "    See 04_diversity/candidate_frame_hotspot_ranking.csv",
  "",
  "Q8. Where are VMEs recorded and which VME types occur?",
  "    See 08_vme_occurrence/vme_presence_across_all_point_frames.csv",
  "    See 08_vme_occurrence/vme_type_occurrence_summary.csv",
  "    See 08_vme_occurrence/vme_presence_by_dive.*",
  "    See 08_vme_occurrence/vme_presence_vs_depth.*",
  "    See 08_vme_occurrence/vme_presence_vs_temperature.*",
  "",
  "CAUTIONS",
  "--------",
  "",
  paste0(
    "- Current community units are: ",
    community_unit_mode,
    "."
  ),
  "- Common-ID richness is standardised annotation-unit richness, not automatically species richness.",
  "- Exclude and UNSURE points are removed from valid ecological denominators.",
  "- VME annotations are analysed separately and never enter point denominators, point relative abundance, richness, or Bray-Curtis matrices.",
  "- VME absence is inferred only across frames represented in point_annotations.csv; VME frames outside that frame universe are reported separately.",
  "- Points with conflicting top-level labels are flagged and excluded.",
  "- Frames may be temporally autocorrelated within dives.",
  "- Environmental models are exploratory and should not be treated as confirmatory.",
  "- Depth and temperature may be correlated; inspect the correlation output.",
  "- Raw dive-level gamma richness is affected by unequal numbers of frames."
)

if (nrow(model_registry_table) > 0) {
  summary_lines <- c(
    summary_lines,
    "",
    "EXPLORATORY GAM RESULTS",
    "-----------------------"
  )

  for (
    row_index in seq_len(
      nrow(model_registry_table)
    )
  ) {
    model_row <-
      model_registry_table[
        row_index,
        ,
        drop = FALSE
      ]

    summary_lines <- c(
      summary_lines,
      paste0(
        "- ",
        model_row$model,
        ": n = ",
        model_row$n_frames,
        ", smooth p = ",
        format(
          model_row$smooth_p_value,
          digits = 4
        ),
        ". Inspect the model plot and diagnostics before interpretation."
      )
    )
  }
}

if (nrow(strong_group_associations) > 0) {
  summary_lines <- c(
    summary_lines,
    "",
    "STRONGEST EXPLORATORY GROUP-ENVIRONMENT CORRELATIONS",
    "---------------------------------------------------"
  )

  for (
    row_index in seq_len(
      nrow(strong_group_associations)
    )
  ) {
    association <-
      strong_group_associations[
        row_index,
        ,
        drop = FALSE
      ]

    summary_lines <- c(
      summary_lines,
      paste0(
        "- ",
        association$community_unit,
        " vs ",
        association$predictor,
        ": Spearman rho = ",
        round(
          association$rho,
          3
        ),
        ", BH-adjusted p = ",
        format(
          association$p_adjusted_bh,
          digits = 4
        ),
        ", n = ",
        association$n,
        "."
      )
    )
  }
}

if (nrow(top_hotspots) > 0) {
  summary_lines <- c(
    summary_lines,
    "",
    "TOP CANDIDATE HOTSPOT FRAMES",
    "------------------------------"
  )

  for (
    row_index in seq_len(
      nrow(top_hotspots)
    )
  ) {
    hotspot <-
      top_hotspots[
        row_index,
        ,
        drop = FALSE
      ]

    summary_lines <- c(
      summary_lines,
      paste0(
        "- ",
        hotspot$dive_id_analysis,
        " / ",
        hotspot$frame_id_analysis,
        ": hotspot score = ",
        round(
          hotspot$hotspot_score,
          3
        ),
        ", rarefied richness = ",
        round(
          hotspot$rarefied_common_id_richness,
          2
        ),
        ", Shannon = ",
        round(
          hotspot$shannon_common_id_diversity,
          3
        ),
        ", biotic cover = ",
        round(
          hotspot$pct_biotic_valid,
          1
        ),
        "%."
      )
    )
  }
}

write_text(
  summary_lines,
  file.path(
    directories$summaries,
    "ecological_questions_and_results.txt"
  )
)


# =============================================================================
# 21. PIPELINE LOG
# =============================================================================

pipeline_log <- c(
  "BIIGLE ecological pipeline completed.",
  "",
  paste0(
    "Point input: ",
    normalizePath(point_input_csv)
  ),
  paste0(
    "VME input: ",
    normalizePath(vme_input_csv)
  ),
  paste0(
    "Output: ",
    normalizePath(output_root)
  ),
  paste0(
    "Point rows read: ",
    format(nrow(point_raw_data), big.mark = ",")
  ),
  paste0(
    "VME rows read: ",
    format(nrow(vme_raw_data), big.mark = ",")
  ),
  paste0(
    "Unique points: ",
    format(
      nrow(point_data),
      big.mark = ","
    )
  ),
  paste0(
    "Point frames: ",
    n_distinct(point_data$sample_id_analysis)
  ),
  paste0(
    "Point frames not exactly ",
    TARGET_POINTS_PER_FRAME,
    " points: ",
    nrow(frames_not_target_point_count)
  ),
  paste0(
    "Point frames with VME: ",
    sum(frame_summary$has_vme)
  ),
  paste0(
    "Unique VME annotations: ",
    nrow(vme_annotations_unique)
  ),
  paste0(
    "VME types: ",
    n_distinct(vme_frame_type_occurrence$vme_unit)
  ),
  paste0(
    "Community-unit mode: ",
    community_unit_mode
  ),
  paste0(
    "Rows with common_id_full: ",
    sum(!is.na(common_id_full_raw))
  ),
  paste0(
    "Rows with depth_m: ",
    sum(!is.na(safe_numeric(column_or_na(raw_data, DEPTH_COL, numeric = TRUE))))
  ),
  paste0(
    "Rows with temperature: ",
    sum(!is.na(safe_numeric(column_or_na(raw_data, TEMPERATURE_COL, numeric = TRUE))))
  ),
  paste0(
    "Frames with depth: ",
    sum(
      !is.na(frame_summary$depth_m)
    )
  ),
  paste0(
    "Frames with temperature: ",
    sum(
      !is.na(frame_summary$temperature)
    )
  ),
  paste0(
    "NMDS stress: ",
    round(
      nmds_stress,
      4
    )
  ),
  paste0(
    "Exploratory GAMs produced: ",
    nrow(model_registry_table)
  )
)

write_text(
  pipeline_log,
  file.path(
    directories$logs,
    "pipeline_summary.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    directories$logs,
    "sessionInfo.txt"
  )
)

message("")
message("BIIGLE ecological analysis complete")
message(
  paste0(
    "Community-unit mode: ",
    community_unit_mode
  )
)
message(
  paste0(
    "Outputs written to: ",
    normalizePath(output_root)
  )
)
message("")
message("Start with:")
message(
  file.path(
    directories$summaries,
    "ecological_questions_and_results.txt"
  )
)
message(
  file.path(
    directories$diversity,
    "frame_alpha_diversity.csv"
  )
)
message(
  file.path(
    directories$ordination,
    "frame_nmds_by_dive.pdf"
  )
)
message(
  file.path(
    directories$vme,
    "vme_presence_across_all_point_frames.csv"
  )
)
