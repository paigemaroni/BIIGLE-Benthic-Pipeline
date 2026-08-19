# =============================================================================
# Shared R helpers for the BIIGLE benthic analysis modules
# =============================================================================

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr", "vegan", "mgcv"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(vegan)
  library(mgcv)
})

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

first_non_missing <- function(x) {
  x <- clean_chr(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_character_ else x[[1]]
}

median_or_na <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

collapse_unique <- function(x, sep = " | ") {
  x <- unique(clean_chr(x))
  x <- sort(x[!is.na(x)])
  if (length(x) == 0) NA_character_ else paste(x, collapse = sep)
}

resolve_single <- function(x) {
  x <- unique(clean_chr(x))
  x <- x[!is.na(x)]
  if (length(x) == 1) x[[1]] else NA_character_
}

n_unique_nonmissing <- function(x) {
  x <- unique(clean_chr(x))
  sum(!is.na(x))
}

standardise_top_level <- function(x) {
  x <- clean_chr(x)
  out <- case_when(
    is.na(x) ~ NA_character_,
    str_to_lower(x) == "biotic" ~ "Biotic",
    str_to_lower(x) == "abiotic" ~ "Abiotic",
    str_to_lower(x) == "exclude" ~ "Exclude",
    str_to_lower(x) == "unsure" ~ "Unsure",
    str_to_lower(x) == "vme" ~ "VME",
    TRUE ~ x
  )
  out
}

resolve_point_class <- function(x) {
  vals <- sort(unique(standardise_top_level(x)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return("Other")
  if (length(vals) == 1) return(vals[[1]])
  return("Conflict")
}

check_required_columns <- function(data, columns, context) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      paste0(context, " is missing required columns: ", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_analysis_dirs <- function() {
  dirs <- c(
    ANALYSES_DIR,
    FIGURES_DIR,
    file.path(ANALYSES_DIR, "00_Quality_Control"),
    file.path(ANALYSES_DIR, "01_Frame_Composition"),
    file.path(ANALYSES_DIR, "02_Alpha_Diversity", "Richness"),
    file.path(ANALYSES_DIR, "02_Alpha_Diversity", "Shannon_Wiener"),
    file.path(ANALYSES_DIR, "02_Alpha_Diversity", "Simpson"),
    file.path(ANALYSES_DIR, "02_Alpha_Diversity", "Evenness"),
    file.path(ANALYSES_DIR, "02_Alpha_Diversity", "Rarefaction"),
    file.path(ANALYSES_DIR, "03_Community_Composition", "Bray_Curtis"),
    file.path(ANALYSES_DIR, "03_Community_Composition", "NMDS"),
    file.path(ANALYSES_DIR, "04_PERMANOVA"),
    file.path(ANALYSES_DIR, "05_CAP_dbRDA"),
    file.path(ANALYSES_DIR, "06_Environmental_Associations"),
    file.path(ANALYSES_DIR, "07_Univariate_Tests"),
    file.path(ANALYSES_DIR, "08_Taxon_Responses"),
    file.path(ANALYSES_DIR, "09_VME")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

save_figure <- function(plot, filename_stub, width = FIGURE_WIDTH, height = FIGURE_HEIGHT) {
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    file.path(FIGURES_DIR, paste0(filename_stub, ".png")),
    plot = plot, width = width, height = height, dpi = FIGURE_DPI
  )
  ggsave(
    file.path(FIGURES_DIR, paste0(filename_stub, ".pdf")),
    plot = plot, width = width, height = height
  )
}

apply_latitude_colour_scale <- function(plot) {
  if (COLOUR_MODE == "custom") {
    if (is.na(LATITUDE_SOUTH_COLOUR) || is.na(LATITUDE_NORTH_COLOUR)) {
      stop(
        "COLOUR_MODE='custom' requires LATITUDE_SOUTH_COLOUR and LATITUDE_NORTH_COLOUR.",
        call. = FALSE
      )
    }
    plot + scale_colour_gradient(
      low = LATITUDE_SOUTH_COLOUR,
      high = LATITUDE_NORTH_COLOUR,
      na.value = "grey70"
    )
  } else {
    plot + scale_colour_viridis_c(
      option = LATITUDE_PALETTE,
      direction = LATITUDE_DIRECTION,
      na.value = "grey70"
    )
  }
}

read_final_point_data <- function() {
  if (!file.exists(POINT_DATA_FILE)) {
    stop(
      paste0(
        "Missing final point dataset: ", POINT_DATA_FILE,
        "\nRun ./Scripts/run_processing_pipeline.sh . first."
      ),
      call. = FALSE
    )
  }
  read_csv(
    POINT_DATA_FILE,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "N/A", "NULL")
  )
}

read_final_vme_data <- function() {
  if (!file.exists(VME_DATA_FILE)) {
    stop(
      paste0(
        "Missing final VME dataset: ", VME_DATA_FILE,
        "\nRun ./Scripts/run_processing_pipeline.sh . first."
      ),
      call. = FALSE
    )
  }
  read_csv(
    VME_DATA_FILE,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "N/A", "NULL")
  )
}

resolve_point_annotations <- function(point_raw) {
  required <- c(
    "dive_id", "filename", "annotation_id", "shape_name", "top_level",
    "label_name", "common_id_short", "common_id_mid", "common_id_full",
    "lat", "long", "depth_m", "temperature",
    "frame_substrate_class", "frame_relief_class"
  )
  check_required_columns(point_raw, required, "point_annotations.csv")

  # The final point file must contain no VME rows.
  if (any(standardise_top_level(point_raw$top_level) == "VME", na.rm = TRUE)) {
    stop("point_annotations.csv contains VME rows.", call. = FALSE)
  }

  point_rows <- point_raw %>%
    mutate(
      .row_id = row_number(),
      shape_std = str_to_lower(clean_chr(shape_name)),
      top_std = standardise_top_level(top_level),
      annotation_key = clean_chr(annotation_id),
      annotation_key = if_else(
        is.na(annotation_key),
        paste0("__missing_annotation_id_row_", .row_id),
        annotation_key
      ),
      community_unit_row = case_when(
        top_std == "Biotic" ~ coalesce(
          clean_chr(common_id_full),
          clean_chr(common_id_mid),
          clean_chr(common_id_short),
          clean_chr(label_name)
        ),
        TRUE ~ NA_character_
      ),
      community_unit_source_row = case_when(
        top_std == "Biotic" & !is.na(clean_chr(common_id_full)) ~ "common_id_full",
        top_std == "Biotic" & !is.na(clean_chr(common_id_mid)) ~ "common_id_mid",
        top_std == "Biotic" & !is.na(clean_chr(common_id_short)) ~ "common_id_short",
        top_std == "Biotic" & !is.na(clean_chr(label_name)) ~ "label_name_fallback",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(shape_std == "point")

  resolved <- point_rows %>%
    group_by(filename, annotation_key) %>%
    summarise(
      annotation_id = first_non_missing(annotation_id),
      dive_id = first_non_missing(dive_id),
      point_class = resolve_point_class(top_std),
      point_class_labels = collapse_unique(top_std),
      community_unit = resolve_single(community_unit_row),
      community_unit_candidates = collapse_unique(community_unit_row),
      community_unit_source = resolve_single(community_unit_source_row),
      multiple_community_units = n_unique_nonmissing(community_unit_row) > 1,
      lat = median_or_na(lat),
      long = median_or_na(long),
      depth_m_raw = median_or_na(depth_m),
      depth_m = if_else(is.na(depth_m_raw), NA_real_, abs(depth_m_raw)),
      temperature = median_or_na(temperature),
      frame_substrate_class = first_non_missing(frame_substrate_class),
      frame_relief_class = first_non_missing(frame_relief_class),
      .groups = "drop"
    )

  resolved
}

build_frame_summary <- function(resolved_points) {
  resolved_points %>%
    group_by(filename, dive_id) %>%
    summarise(
      n_points = n(),
      n_biotic = sum(point_class == "Biotic", na.rm = TRUE),
      n_abiotic = sum(point_class == "Abiotic", na.rm = TRUE),
      n_exclude = sum(point_class == "Exclude", na.rm = TRUE),
      n_unsure = sum(point_class == "Unsure", na.rm = TRUE),
      n_conflict = sum(point_class == "Conflict", na.rm = TRUE),
      n_other = sum(!point_class %in% c(
        "Biotic", "Abiotic", "Exclude", "Unsure", "Conflict"
      ), na.rm = TRUE),
      n_biotic_assigned = sum(
        point_class == "Biotic" &
          !is.na(community_unit) &
          !multiple_community_units,
        na.rm = TRUE
      ),
      n_biotic_unassigned = sum(
        point_class == "Biotic" &
          (is.na(community_unit) | multiple_community_units),
        na.rm = TRUE
      ),
      lat = median_or_na(lat),
      long = median_or_na(long),
      depth_m = median_or_na(depth_m),
      depth_m_raw = median_or_na(depth_m_raw),
      temperature = median_or_na(temperature),
      frame_substrate_class = first_non_missing(frame_substrate_class),
      frame_relief_class = first_non_missing(frame_relief_class),
      .groups = "drop"
    ) %>%
    mutate(
      target_points = TARGET_POINTS_PER_FRAME,
      difference_from_target = n_points - target_points,
      point_count_qc = if_else(n_points == target_points, "PASS", "REVIEW"),
      n_valid = n_biotic + n_abiotic,
      pct_biotic_all = if_else(n_points > 0, 100 * n_biotic / n_points, NA_real_),
      pct_abiotic_all = if_else(n_points > 0, 100 * n_abiotic / n_points, NA_real_),
      pct_exclude_all = if_else(n_points > 0, 100 * n_exclude / n_points, NA_real_),
      pct_unsure_all = if_else(n_points > 0, 100 * n_unsure / n_points, NA_real_),
      pct_biotic_valid = if_else(n_valid > 0, 100 * n_biotic / n_valid, NA_real_),
      pct_abiotic_valid = if_else(n_valid > 0, 100 * n_abiotic / n_valid, NA_real_),
      pct_biotic_units_resolved = if_else(
        n_biotic > 0, 100 * n_biotic_assigned / n_biotic, NA_real_
      )
    )
}

build_community_counts <- function(resolved_points) {
  resolved_points %>%
    filter(
      point_class == "Biotic",
      !is.na(community_unit),
      !multiple_community_units
    ) %>%
    count(filename, dive_id, community_unit, name = "abundance")
}

community_matrix_from_counts <- function(count_table, row_id = "filename") {
  if (nrow(count_table) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  formula <- as.formula(paste0("abundance ~ ", row_id, " + community_unit"))
  mat <- xtabs(formula, data = count_table)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

transform_community_matrix <- function(mat, method = COMMUNITY_TRANSFORM) {
  method <- tolower(method)

  if (method == "none") return(mat)
  if (method == "sqrt") return(sqrt(mat))
  if (method == "relative") return(decostand(mat, method = "total"))
  if (method == "sqrt_relative") {
    return(sqrt(decostand(mat, method = "total")))
  }

  stop(
    paste0(
      "Unknown COMMUNITY_TRANSFORM: ", method,
      ". Use none, sqrt, relative, or sqrt_relative."
    ),
    call. = FALSE
  )
}



community_model_keep <- function(frame_summary, community_matrix) {
  frame_summary %>%
    filter(
      filename %in% rownames(community_matrix),
      n_points >= MODEL_MIN_TOTAL_POINTS,
      n_points <= MODEL_MAX_TOTAL_POINTS,
      n_biotic_assigned >= MIN_BIOTIC_POINTS
    ) %>%
    pull(filename)
}

frame_alpha_diversity <- function(community_matrix, frame_summary) {
  # Alpha diversity is reported against the COMPLETE point-sampled frame
  # universe, including frames with zero resolved Biotic points.
  #
  # For zero-Biotic frames:
  #   observed richness = 0 (nothing biological was observed at random points)
  #   Shannon / inverse Simpson / Pielou = NA (no assemblage distribution exists)
  #   rarefied richness = NA (insufficient biological observations)
  #
  # Inferential alpha-diversity models should use alpha_model_eligible, which
  # enforces both the configured overall point-effort range and minimum number
  # of resolved Biotic points.

  base <- frame_summary %>%
    mutate(
      assigned_biotic_points = n_biotic_assigned,
      observed_common_id_richness = 0,
      rarefied_common_id_richness = NA_real_,
      shannon_common_id_diversity = NA_real_,
      inverse_simpson_common_id_diversity = NA_real_,
      pielou_evenness = NA_real_
    )

  if (nrow(community_matrix) > 0) {
    totals <- rowSums(community_matrix)
    richness <- rowSums(community_matrix > 0)
    shannon <- diversity(community_matrix, index = "shannon")
    invsimpson <- diversity(community_matrix, index = "invsimpson")
    evenness <- ifelse(richness > 1, shannon / log(richness), NA_real_)

    rarefied <- rep(NA_real_, nrow(community_matrix))
    eligible <- totals >= RAREFACTION_N

    if (any(eligible)) {
      rarefied[eligible] <- rarefy(
        community_matrix[eligible, , drop = FALSE],
        sample = RAREFACTION_N
      )
    }

    calculated <- tibble(
      filename = rownames(community_matrix),
      assigned_biotic_points_calc = totals,
      observed_common_id_richness_calc = richness,
      rarefied_common_id_richness_calc = rarefied,
      shannon_common_id_diversity_calc = shannon,
      inverse_simpson_common_id_diversity_calc = invsimpson,
      pielou_evenness_calc = evenness
    )

    base <- base %>%
      left_join(calculated, by = "filename") %>%
      mutate(
        assigned_biotic_points = coalesce(
          assigned_biotic_points_calc,
          assigned_biotic_points
        ),
        observed_common_id_richness = coalesce(
          observed_common_id_richness_calc,
          observed_common_id_richness
        ),
        rarefied_common_id_richness = rarefied_common_id_richness_calc,
        shannon_common_id_diversity = shannon_common_id_diversity_calc,
        inverse_simpson_common_id_diversity = inverse_simpson_common_id_diversity_calc,
        pielou_evenness = pielou_evenness_calc
      ) %>%
      select(-ends_with("_calc"))
  }

  base %>%
    mutate(
      has_resolved_biota = assigned_biotic_points > 0,
      rarefaction_eligible = assigned_biotic_points >= RAREFACTION_N,
      alpha_model_eligible =
        n_points >= MODEL_MIN_TOTAL_POINTS &
        n_points <= MODEL_MAX_TOTAL_POINTS &
        assigned_biotic_points >= MIN_BIOTIC_POINTS,
      alpha_model_exclusion_reason = case_when(
        n_points < MODEL_MIN_TOTAL_POINTS ~ paste0(
          "total_points_below_", MODEL_MIN_TOTAL_POINTS
        ),
        n_points > MODEL_MAX_TOTAL_POINTS ~ paste0(
          "total_points_above_", MODEL_MAX_TOTAL_POINTS
        ),
        assigned_biotic_points < MIN_BIOTIC_POINTS ~ paste0(
          "biotic_points_below_", MIN_BIOTIC_POINTS
        ),
        TRUE ~ "eligible"
      )
    )
}

prepare_point_analysis <- function() {
  ensure_analysis_dirs()
  raw <- read_final_point_data()
  resolved <- resolve_point_annotations(raw)
  frame_summary <- build_frame_summary(resolved)
  counts <- build_community_counts(resolved)
  community <- community_matrix_from_counts(counts, "filename")
  transformed <- transform_community_matrix(community)

  list(
    raw = raw,
    resolved_points = resolved,
    frame_summary = frame_summary,
    community_counts = counts,
    community_matrix = community,
    transformed_community_matrix = transformed
  )
}

write_matrix_csv <- function(mat, path, row_name = "filename") {
  if (nrow(mat) == 0) {
    write_csv(tibble(), path)
    return(invisible(NULL))
  }
  out <- as.data.frame(mat, check.names = FALSE)
  out[[row_name]] <- rownames(mat)
  out <- out[, c(row_name, setdiff(names(out), row_name)), drop = FALSE]
  write_csv(out, path)
}

extract_anova_table <- function(x) {
  out <- as.data.frame(x)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("term", setdiff(names(out), "term")), drop = FALSE]
  as_tibble(out)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(safe_num(x)) & is.finite(safe_num(y))
  if (sum(keep) < 4) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
  }
  test <- suppressWarnings(cor.test(
    safe_num(x)[keep], safe_num(y)[keep],
    method = "spearman", exact = FALSE
  ))
  tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n = sum(keep)
  )
}

model_sample_report <- function(data, variables, analysis_name) {
  present <- intersect(variables, names(data))
  complete <- if (length(present) == 0) rep(TRUE, nrow(data)) else complete.cases(data[, present, drop = FALSE])
  tibble(
    analysis = analysis_name,
    total_rows_available = nrow(data),
    rows_complete_for_model = sum(complete),
    rows_excluded_missing = sum(!complete),
    variables = paste(present, collapse = " + ")
  )
}
