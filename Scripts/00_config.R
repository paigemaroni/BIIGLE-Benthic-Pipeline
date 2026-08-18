# =============================================================================
# BIIGLE BENTHIC PIPELINE — CENTRAL CONFIGURATION
# =============================================================================
# All user-editable analysis settings should live here.
# Analysis scripts source this file so visual/statistical changes can be made
# without editing the analysis logic itself.

PROJECT_ROOT <- "."
SHEETS_DIR   <- file.path(PROJECT_ROOT, "Sheets")
DATA_DIR     <- file.path(PROJECT_ROOT, "Data")
ANALYSES_DIR <- file.path(PROJECT_ROOT, "Analyses")
FIGURES_DIR  <- file.path(PROJECT_ROOT, "Figures")

POINT_DATA_FILE <- file.path(DATA_DIR, "Final", "point_annotations.csv")
VME_DATA_FILE   <- file.path(DATA_DIR, "Final", "vme_annotations.csv")

# Sampling / diversity
TARGET_POINTS_PER_FRAME <- 100
MIN_BIOTIC_POINTS <- 10
MODEL_MIN_TOTAL_POINTS <- 80
MODEL_MAX_TOTAL_POINTS <- 120
RAREFACTION_N <- 10
TOP_N_TAXA <- 12

# Community resemblance / ordination
COMMUNITY_TRANSFORM <- "sqrt_relative"
# Allowed: "none", "sqrt", "relative", "sqrt_relative"
DISTANCE_METHOD <- "bray"
NMDS_DIMENSIONS <- 2
NMDS_TRYMAX <- 200
N_PERMUTATIONS <- 9999
PERMUTATION_STRATA_COLUMN <- "dive_id"

# Candidate environmental variables.
CONTINUOUS_ENVIRONMENTAL_VARIABLES <- c(
  "depth_m", "temperature", "lat", "long"
)

CATEGORICAL_ENVIRONMENTAL_VARIABLES <- c(
  "frame_substrate_class", "frame_relief_class"
)

ENVIRONMENTAL_VARIABLES <- c(
  CONTINUOUS_ENVIRONMENTAL_VARIABLES,
  CATEGORICAL_ENVIRONMENTAL_VARIABLES
)

# Frame-level PERMANOVA terms. These are tested with permutations restricted
# within dive where possible. Latitude/longitude are handled at dive level
# because they are usually constant within a dive.
PERMANOVA_FRAME_TERMS <- c(
  "depth_m", "temperature",
  "frame_substrate_class", "frame_relief_class"
)

PERMANOVA_DIVE_TERMS <- c("lat", "long")

# Univariate screening
MIN_GROUP_N <- 5
P_ADJUST_METHOD <- "BH"
ALPHA <- 0.05

# Plot settings
FIGURE_DPI <- 300
FIGURE_WIDTH <- 9
FIGURE_HEIGHT <- 7
POINT_SIZE <- 2.7
POINT_ALPHA <- 0.80

# Continuous latitude palette:
# higher/northern latitude should render darker by default.
COLOUR_MODE <- "palette"       # "palette" or "custom"
LATITUDE_PALETTE <- "viridis"  # ggplot2 viridis option name
LATITUDE_DIRECTION <- -1
LATITUDE_SOUTH_COLOUR <- NA_character_
LATITUDE_NORTH_COLOUR <- NA_character_

SHOW_NMDS_LABELS <- FALSE
SHOW_NMDS_ELLIPSES <- TRUE
NMDS_ELLIPSE_GROUP <- "frame_substrate_class"

# Reproducibility
RANDOM_SEED <- 20260818
set.seed(RANDOM_SEED)
