#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 22: Bray-Curtis and NMDS\n")
cat("=================================\n")

a <- prepare_point_analysis()
mat <- a$community_matrix
meta <- a$frame_summary

eligible_frames <- community_model_keep(meta, mat)
mat <- mat[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]

if (nrow(mat) < 4 || ncol(mat) < 2) {
  stop("Too few eligible frames or biological groups for NMDS.", call. = FALSE)
}

mat_t <- transform_community_matrix(mat, COMMUNITY_TRANSFORM)
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

set.seed(RANDOM_SEED)
nmds <- metaMDS(
  mat_t,
  distance = DISTANCE_METHOD,
  k = NMDS_DIMENSIONS,
  trymax = NMDS_TRYMAX,
  autotransform = FALSE,
  trace = FALSE
)

site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$filename <- rownames(site_scores)

scores_out <- as_tibble(site_scores) %>%
  left_join(meta, by = "filename")

write_csv(scores_out, file.path(nmds_dir, "22_nmds_scores.csv"))

write_csv(
  tibble(
    metric = c(
      "stress", "dimensions", "frames", "biological_groups",
      "distance_method", "community_transform", "minimum_biotic_points"
    ),
    value = c(
      format(nmds$stress, digits = 8),
      NMDS_DIMENSIONS,
      nrow(mat_t),
      ncol(mat_t),
      DISTANCE_METHOD,
      COMMUNITY_TRANSFORM,
      MIN_BIOTIC_POINTS
    )
  ),
  file.path(nmds_dir, "22_nmds_model_summary.csv")
)

axis_x <- names(site_scores)[1]
axis_y <- names(site_scores)[2]

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

  if (
    SHOW_NMDS_ELLIPSES &&
    variable == NMDS_ELLIPSE_GROUP
  ) {
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

# Fit environmental variables separately so missing depth/temperature does not
# discard frames from latitude/longitude fits.
envfit_rows <- list()

for (variable in ENVIRONMENTAL_VARIABLES) {
  if (!variable %in% names(scores_out)) next

  vals <- scores_out[[variable]]
  keep_var <- !is.na(vals) & as.character(vals) != ""

  if (sum(keep_var) < 4) next

  score_matrix <- as.matrix(scores_out[keep_var, c(axis_x, axis_y), drop = FALSE])
  env_data <- data.frame(value = vals[keep_var])

  if (variable %in% CONTINUOUS_ENVIRONMENTAL_VARIABLES) {
    env_data$value <- safe_num(env_data$value)
    good <- is.finite(env_data$value)
    score_matrix <- score_matrix[good, , drop = FALSE]
    env_data <- env_data[good, , drop = FALSE]
    if (nrow(env_data) < 4 || length(unique(env_data$value)) < 3) next
  } else {
    env_data$value <- factor(env_data$value)
    if (nlevels(env_data$value) < 2) next
  }

  set.seed(RANDOM_SEED)
  fit <- envfit(
    score_matrix,
    env_data,
    permutations = N_PERMUTATIONS,
    na.rm = TRUE
  )

  if (variable %in% CONTINUOUS_ENVIRONMENTAL_VARIABLES) {
    vec <- scores(fit, display = "vectors")
    envfit_rows[[length(envfit_rows) + 1]] <- tibble(
      variable = variable,
      variable_type = "continuous",
      axis1 = vec[1, 1],
      axis2 = vec[1, 2],
      r_squared = unname(fit$vectors$r[[1]]),
      p_value = unname(fit$vectors$pvals[[1]]),
      n = nrow(env_data)
    )
  } else {
    envfit_rows[[length(envfit_rows) + 1]] <- tibble(
      variable = variable,
      variable_type = "factor",
      axis1 = NA_real_,
      axis2 = NA_real_,
      r_squared = unname(fit$factors$r[[1]]),
      p_value = unname(fit$factors$pvals[[1]]),
      n = nrow(env_data)
    )
  }
}

if (length(envfit_rows) > 0) {
  envfit_table <- bind_rows(envfit_rows) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(envfit_table, file.path(nmds_dir, "22_nmds_environmental_fit.csv"))
}

saveRDS(nmds, file.path(nmds_dir, "22_nmds_model.rds"))

cat("NMDS frames: ", nrow(mat_t), "\n", sep = "")
cat("NMDS groups: ", ncol(mat_t), "\n", sep = "")
cat("Stress: ", format(nmds$stress, digits = 5), "\n", sep = "")
