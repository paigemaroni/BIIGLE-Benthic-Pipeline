# =============================================================================
# BIIGLE BENTHIC PIPELINE — CENTRAL CONFIGURATION
# =============================================================================
# Future refactored R analysis modules should source this file.

PROJECT_ROOT <- "."
SHEETS_DIR   <- file.path(PROJECT_ROOT, "Sheets")
DATA_DIR     <- file.path(PROJECT_ROOT, "Data")
ANALYSES_DIR <- file.path(PROJECT_ROOT, "Analyses")
FIGURES_DIR  <- file.path(PROJECT_ROOT, "Figures")

TARGET_POINTS_PER_FRAME <- 100
MIN_BIOTIC_POINTS <- 10
RAREFACTION_N <- 10

COMMUNITY_TRANSFORM <- "sqrt"
DISTANCE_METHOD <- "bray"
NMDS_DIMENSIONS <- 2
NMDS_TRYMAX <- 200
N_PERMUTATIONS <- 9999
PERMUTATION_STRATA_COLUMN <- "dive_id"

ENVIRONMENTAL_VARIABLES <- c(
  "depth_m", "temperature", "lat", "long",
  "frame_substrate_class", "frame_relief_class"
)

P_ADJUST_METHOD <- "BH"
ALPHA <- 0.05

# Plot settings
FIGURE_DPI <- 300
FIGURE_WIDTH <- 9
FIGURE_HEIGHT <- 7
POINT_SIZE <- 2.7
POINT_ALPHA <- 0.80

# Continuous latitude palette: higher latitude should render darker by default.
COLOUR_MODE <- "palette"       # "palette" or "custom"
LATITUDE_PALETTE <- "viridis"
LATITUDE_DIRECTION <- -1
LATITUDE_SOUTH_COLOUR <- NA_character_
LATITUDE_NORTH_COLOUR <- NA_character_

SHOW_NMDS_LABELS <- FALSE
SHOW_NMDS_ELLIPSES <- TRUE
NMDS_ELLIPSE_GROUP <- "frame_substrate_class"

RANDOM_SEED <- 20260817
set.seed(RANDOM_SEED)
