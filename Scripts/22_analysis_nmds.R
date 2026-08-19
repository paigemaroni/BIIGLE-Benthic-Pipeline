#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 22: Bray-Curtis and NMDS\n")
cat("=================================\n")

a <- prepare_point_analysis()
mat_all <- a$community_matrix
meta_all <- a$frame_summary

eligible_frames <- community_model_keep(meta_all, mat_all)
mat <- mat_all[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]
meta <- meta_all %>% filter(filename %in% rownames(mat))
meta <- meta[match(rownames(mat), meta$filename), , drop = FALSE]

if (nrow(mat) < 4 || ncol(mat) < 2) {
  stop("Too few eligible frames or biological groups for NMDS.", call. = FALSE)
}

mat_t <- force_numeric_matrix(
  transform_community_matrix(mat, COMMUNITY_TRANSFORM),
  "NMDS transformed community matrix"
)

if (any(!is.finite(mat_t))) {
  stop("NMDS transformed community matrix contains non-finite values.", call. = FALSE)
}
if (any(mat_t < 0)) {
  stop("NMDS transformed community matrix contains negative values.", call. = FALSE)
}
if (any(rowSums(mat_t) <= 0)) {
  stop("NMDS transformed community matrix contains one or more empty rows.", call. = FALSE)
}

bray <- vegdist(mat_t, method = DISTANCE_METHOD)

bray_dir <- file.path(ANALYSES_DIR, "03_Community_Composition", "Bray_Curtis")
nmds_dir <- file.path(ANALYSES_DIR, "03_Community_Composition", "NMDS")
dir.create(bray_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nmds_dir, recursive = TRUE, showWarnings = FALSE)

write_matrix_csv(
  as.matrix(bray),
  file.path(bray_dir, "22_bray_curtis_matrix.csv"),
  row_name = "filename"
)

# Bray-Curtis = 1 means complete dissimilarity; for non-negative community
# data this commonly reflects pairs with no shared community units. Report it
# explicitly because many no-share pairs can make 2D NMDS difficult to fit.
bray_values <- as.numeric(bray)
n_pairs <- length(bray_values)
n_complete_dissimilarity <- sum(bray_values >= (1 - 1e-12), na.rm = TRUE)
pct_complete_dissimilarity <- if (n_pairs > 0) {
  100 * n_complete_dissimilarity / n_pairs
} else {
  NA_real_
}

set.seed(RANDOM_SEED)
nmds <- metaMDS(
  bray,
  k = NMDS_DIMENSIONS,
  trymax = NMDS_TRYMAX,
  autotransform = FALSE,
  trace = FALSE
)

# Diagnostic dimensionality check. A weak 2-D solution should not be forced
# into interpretation when a low-dimensional 3-D solution represents the
# same Bray-Curtis ranks substantially better.
nmds_3d <- NULL
if (NMDS_DIMENSIONS == 2 && nrow(mat_t) >= 5 && ncol(mat_t) >= 3) {
  set.seed(RANDOM_SEED)
  nmds_3d <- metaMDS(
    bray,
    k = 3,
    trymax = NMDS_TRYMAX,
    autotransform = FALSE,
    trace = FALSE
  )

  stress_reduction_pct <- if (is.finite(nmds$stress) && nmds$stress > 0) {
    100 * (nmds$stress - nmds_3d$stress) / nmds$stress
  } else {
    NA_real_
  }

  dimension_comparison <- tibble(
    dimensions = c(2L, 3L),
    stress = c(nmds$stress, nmds_3d$stress),
    relative_stress_reduction_from_2d_percent = c(0, stress_reduction_pct)
  )

  write_csv(
    dimension_comparison,
    file.path(nmds_dir, "22_nmds_dimension_comparison.csv")
  )

  site_scores_3d <- as.data.frame(scores(nmds_3d, display = "sites"))
  site_scores_3d$filename <- rownames(site_scores_3d)

  write_csv(
    as_tibble(site_scores_3d) %>% left_join(meta, by = "filename"),
    file.path(nmds_dir, "22_nmds_scores_3d.csv")
  )

  saveRDS(
    nmds_3d,
    file.path(nmds_dir, "22_nmds_model_3d.rds")
  )
}

site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$filename <- rownames(site_scores)

scores_out <- as_tibble(site_scores) %>%
  left_join(meta, by = "filename")

write_csv(scores_out, file.path(nmds_dir, "22_nmds_scores.csv"))

stress_class <- dplyr::case_when(
  nmds$stress < 0.05 ~ "excellent",
  nmds$stress < 0.10 ~ "very_good",
  nmds$stress < 0.20 ~ "potentially_useful_with_caution",
  TRUE ~ "weak_two_dimensional_representation"
)

write_csv(
  tibble(
    metric = c(
      "stress", "stress_interpretation", "dimensions",
      "point_sampled_frames", "eligible_frames", "excluded_frames",
      "biological_groups", "distance_method", "community_transform",
      "minimum_biotic_points", "minimum_total_points", "maximum_total_points",
      "bray_curtis_pairs", "complete_dissimilarity_pairs",
      "complete_dissimilarity_pairs_percent"
    ),
    value = c(
      format(nmds$stress, digits = 8),
      stress_class,
      NMDS_DIMENSIONS,
      nrow(meta_all),
      nrow(mat_t),
      nrow(meta_all) - nrow(mat_t),
      ncol(mat_t),
      DISTANCE_METHOD,
      COMMUNITY_TRANSFORM,
      MIN_BIOTIC_POINTS,
      MODEL_MIN_TOTAL_POINTS,
      MODEL_MAX_TOTAL_POINTS,
      n_pairs,
      n_complete_dissimilarity,
      format(pct_complete_dissimilarity, digits = 8)
    )
  ),
  file.path(nmds_dir, "22_nmds_model_summary.csv")
)

axis_x <- names(site_scores)[1]
axis_y <- names(site_scores)[2]

# Latitude is useful as a visual gradient at frame level, but its significance
# is NOT tested at frame level because latitude is effectively repeated within
# dive. Dive-level inference is handled below using dive centroids.
lat_plot_data <- scores_out %>% filter(!is.na(lat))

if (nrow(lat_plot_data) > 0) {
  p6 <- ggplot(
    lat_plot_data,
    aes(x = .data[[axis_x]], y = .data[[axis_y]], colour = lat)
  ) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    coord_equal() +
    labs(
      title = "Community composition across latitude",
      subtitle = paste0(
        "NMDS of ", COMMUNITY_TRANSFORM, "-transformed abundance; ",
        "Bray-Curtis; stress = ", round(nmds$stress, 3)
      ),
      x = axis_x,
      y = axis_y,
      colour = "Latitude"
    ) +
    theme_bw(base_size = 11)

  p6 <- apply_latitude_colour_scale(p6)
  save_figure(p6, "Fig_06_NMDS_Latitude")
}

plot_factor_nmds <- function(variable, title, file_stub) {
  dat <- scores_out %>% filter(!is.na(.data[[variable]]), .data[[variable]] != "")
  if (nrow(dat) == 0 || n_distinct(dat[[variable]]) < 2) return(invisible(NULL))

  p <- ggplot(
    dat,
    aes(x = .data[[axis_x]], y = .data[[axis_y]], colour = .data[[variable]])
  ) +
    geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
    coord_equal() +
    labs(
      title = title,
      subtitle = paste0("Bray-Curtis NMDS; stress = ", round(nmds$stress, 3)),
      x = axis_x,
      y = axis_y,
      colour = variable
    ) +
    theme_bw(base_size = 10)

  if (SHOW_NMDS_ELLIPSES && variable == NMDS_ELLIPSE_GROUP) {
    group_counts <- dat %>% count(.data[[variable]])
    eligible <- group_counts %>% filter(n >= 3) %>% pull(1)
    ellipse_dat <- dat %>% filter(.data[[variable]] %in% eligible)
    if (n_distinct(ellipse_dat[[variable]]) >= 1) {
      p <- p + stat_ellipse(
        data = ellipse_dat,
        aes(group = .data[[variable]], colour = .data[[variable]]),
        type = "t",
        linewidth = 0.6,
        show.legend = FALSE
      )
    }
  }

  save_figure(p, file_stub, width = 10, height = 8)
}

plot_factor_nmds(
  "frame_substrate_class",
  "Community composition by frame substrate class",
  "Fig_07_NMDS_Substrate"
)

plot_factor_nmds(
  "frame_relief_class",
  "Community composition by frame relief class",
  "Fig_08_NMDS_Relief"
)

# -------------------------------------------------------------------------
# Frame-level environmental fits
# -------------------------------------------------------------------------
# Depth, temperature, substrate and relief may vary among frames. Their
# permutation tests are restricted within dive to respect repeated sampling.
# Latitude and longitude are deliberately excluded here because they are
# effectively dive-level predictors in this worked example.
frame_envfit_variables <- intersect(
  c("depth_m", "temperature", "frame_substrate_class", "frame_relief_class"),
  names(scores_out)
)

frame_envfit_rows <- list()

for (variable in frame_envfit_variables) {
  vals <- scores_out[[variable]]
  keep_var <- !is.na(vals) & as.character(vals) != "" & !is.na(scores_out$dive_id)

  dat <- scores_out[keep_var, , drop = FALSE]
  if (nrow(dat) < 4) next

  score_matrix <- as.matrix(dat[, c(axis_x, axis_y), drop = FALSE])
  env_data <- data.frame(value = dat[[variable]])

  if (variable %in% CONTINUOUS_ENVIRONMENTAL_VARIABLES) {
    env_data$value <- safe_num(env_data$value)
    good <- is.finite(env_data$value)
    score_matrix <- score_matrix[good, , drop = FALSE]
    env_data <- env_data[good, , drop = FALSE]
    dat <- dat[good, , drop = FALSE]
    if (nrow(env_data) < 4 || length(unique(env_data$value)) < 3) next
  } else {
    env_data$value <- droplevels(factor(env_data$value))
    if (nlevels(env_data$value) < 2) next
  }

  set.seed(RANDOM_SEED)
  fit <- envfit(
    score_matrix,
    env_data,
    permutations = N_PERMUTATIONS,
    strata = dat$dive_id,
    na.rm = TRUE
  )

  if (variable %in% CONTINUOUS_ENVIRONMENTAL_VARIABLES) {
    vec <- scores(fit, display = "vectors")
    frame_envfit_rows[[length(frame_envfit_rows) + 1]] <- tibble(
      variable = variable,
      variable_type = "continuous",
      axis1 = vec[1, 1],
      axis2 = vec[1, 2],
      r_squared = unname(fit$vectors$r[[1]]),
      p_value = unname(fit$vectors$pvals[[1]]),
      n = nrow(env_data),
      n_dives = n_distinct(dat$dive_id),
      permutation_restriction = PERMUTATION_STRATA_COLUMN
    )
  } else {
    frame_envfit_rows[[length(frame_envfit_rows) + 1]] <- tibble(
      variable = variable,
      variable_type = "factor",
      axis1 = NA_real_,
      axis2 = NA_real_,
      r_squared = unname(fit$factors$r[[1]]),
      p_value = unname(fit$factors$pvals[[1]]),
      n = nrow(env_data),
      n_dives = n_distinct(dat$dive_id),
      permutation_restriction = PERMUTATION_STRATA_COLUMN
    )
  }
}

if (length(frame_envfit_rows) > 0) {
  frame_envfit_table <- bind_rows(frame_envfit_rows) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))

  write_csv(
    frame_envfit_table,
    file.path(nmds_dir, "22_nmds_environmental_fit_frame_level.csv")
  )
}

# -------------------------------------------------------------------------
# Dive-centroid geographic fits
# -------------------------------------------------------------------------
# Average NMDS coordinates within dive and test latitude/longitude using the
# dive (not the frame) as the inferential unit. These are exploratory
# ordination overlays; formal geographic community tests are performed later
# using dive-level PERMANOVA/dbRDA.
dive_centroids <- scores_out %>%
  filter(!is.na(dive_id)) %>%
  group_by(dive_id) %>%
  summarise(
    axis1 = mean(.data[[axis_x]], na.rm = TRUE),
    axis2 = mean(.data[[axis_y]], na.rm = TRUE),
    lat = median_or_na(lat),
    long = median_or_na(long),
    eligible_frames = n(),
    .groups = "drop"
  )

write_csv(
  dive_centroids,
  file.path(nmds_dir, "22_nmds_dive_centroids.csv")
)

dive_envfit_rows <- list()
for (variable in intersect(c("lat", "long"), names(dive_centroids))) {
  dat <- dive_centroids %>%
    filter(
      is.finite(.data[[variable]]),
      is.finite(axis1),
      is.finite(axis2)
    )

  if (nrow(dat) < 6 || n_distinct(dat[[variable]]) < 3) next

  score_matrix <- as.matrix(dat[, c("axis1", "axis2"), drop = FALSE])
  env_data <- data.frame(value = dat[[variable]])

  set.seed(RANDOM_SEED)
  fit <- envfit(
    score_matrix,
    env_data,
    permutations = N_PERMUTATIONS,
    na.rm = TRUE
  )

  vec <- scores(fit, display = "vectors")
  dive_envfit_rows[[length(dive_envfit_rows) + 1]] <- tibble(
    variable = variable,
    variable_type = "continuous",
    axis1 = vec[1, 1],
    axis2 = vec[1, 2],
    r_squared = unname(fit$vectors$r[[1]]),
    p_value = unname(fit$vectors$pvals[[1]]),
    n_dives = nrow(dat),
    permutation_restriction = "none_dive_is_unit"
  )
}

if (length(dive_envfit_rows) > 0) {
  dive_envfit_table <- bind_rows(dive_envfit_rows) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))

  write_csv(
    dive_envfit_table,
    file.path(nmds_dir, "22_nmds_environmental_fit_dive_level.csv")
  )
}

saveRDS(nmds, file.path(nmds_dir, "22_nmds_model.rds"))

cat("Point-sampled frames: ", nrow(meta_all), "\n", sep = "")
cat("NMDS-eligible frames: ", nrow(mat_t), "\n", sep = "")
cat("Excluded by community-model QC: ", nrow(meta_all) - nrow(mat_t), "\n", sep = "")
cat("NMDS groups: ", ncol(mat_t), "\n", sep = "")
cat("Stress: ", format(nmds$stress, digits = 5), " (", stress_class, ")\n", sep = "")
cat(
  "Bray-Curtis pairs at complete dissimilarity: ",
  n_complete_dissimilarity, "/", n_pairs,
  " (", round(pct_complete_dissimilarity, 2), "%)\n",
  sep = ""
)
if (!is.null(nmds_3d)) {
  cat(
    "3-D diagnostic stress: ",
    format(nmds_3d$stress, digits = 5),
    " (see 22_nmds_dimension_comparison.csv)\n",
    sep = ""
  )
}
cat("Frame-level envfit permutations: restricted within dive\n")
cat("Latitude/longitude envfit unit: dive centroid\n")
