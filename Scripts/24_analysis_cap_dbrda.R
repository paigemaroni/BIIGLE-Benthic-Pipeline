#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 24: CAP / distance-based RDA\n")
cat("=====================================\n")

a <- prepare_point_analysis()
mat_all <- a$community_matrix
meta_all <- a$frame_summary

eligible_frames <- community_model_keep(meta_all, mat_all)
mat <- mat_all[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]
meta <- meta_all %>% filter(filename %in% rownames(mat))
meta <- meta[match(rownames(mat), meta$filename), , drop = FALSE]

out_dir <- file.path(ANALYSES_DIR, "05_CAP_dbRDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trim_rare_factor_levels <- function(dat, terms, min_n = MIN_GROUP_N) {
  factor_terms <- intersect(terms, CATEGORICAL_ENVIRONMENTAL_VARIABLES)
  if (length(factor_terms) == 0 || nrow(dat) == 0) return(dat)

  repeat {
    n_before <- nrow(dat)
    for (variable in factor_terms) {
      vals <- clean_chr(dat[[variable]])
      counts <- table(vals, useNA = "no")
      keep_levels <- names(counts[counts >= min_n])

      dat <- dat[
        !is.na(vals) & vals %in% keep_levels,
        ,
        drop = FALSE
      ]
      if (nrow(dat) == 0) return(dat)
      dat[[variable]] <- droplevels(factor(dat[[variable]]))
    }
    if (nrow(dat) == n_before) break
  }

  dat
}

add_adjusted_p <- function(tab) {
  out <- extract_anova_table(tab)
  if ("Pr(>F)" %in% names(out)) {
    out <- out %>%
      mutate(
        p_value = .data[["Pr(>F)"]],
        p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD)
      )
  }
  out
}

# -------------------------------------------------------------------------
# Frame-level constrained ordination
# -------------------------------------------------------------------------
terms <- intersect(
  c("depth_m", "temperature", "frame_substrate_class", "frame_relief_class"),
  names(meta)
)

if (length(terms) == 0) {
  stop("No configured environmental variables are present for CAP/dbRDA.", call. = FALSE)
}

complete <- complete.cases(meta[, terms, drop = FALSE])
dat <- meta[complete, , drop = FALSE]
n_complete_before_replication_filter <- nrow(dat)

dat <- trim_rare_factor_levels(dat, terms, MIN_GROUP_N)

for (variable in intersect(terms, CATEGORICAL_ENVIRONMENTAL_VARIABLES)) {
  if (variable %in% names(dat)) {
    dat[[variable]] <- droplevels(factor(dat[[variable]]))
  }
}

usable_terms <- terms[
  vapply(
    terms,
    function(variable) {
      variable %in% names(dat) &&
        nrow(dat) > 0 &&
        n_distinct(dat[[variable]]) >= 2
    },
    logical(1)
  )
]

if (nrow(dat) < 15 || length(usable_terms) == 0) {
  stop(
    "Too few replicated complete frames for the configured CAP/dbRDA model.",
    call. = FALSE
  )
}

comm <- force_numeric_matrix(
  transform_community_matrix(
    mat[dat$filename, , drop = FALSE],
    COMMUNITY_TRANSFORM
  ),
  "dbRDA transformed community matrix"
)

formula <- as.formula(
  paste("comm ~", paste(usable_terms, collapse = " + "))
)

# vegan's dbrda() provides distance-based RDA. Permutations are blocked by
# dive using a permute design so repeated frames never exchange among dives.
fit <- dbrda(
  formula,
  data = dat,
  distance = DISTANCE_METHOD
)


# -------------------------------------------------------------------------
# Model diagnostics: explanatory power, collinearity and aliasing
# -------------------------------------------------------------------------
# RsquareAdj() is supported for dbRDA objects. vif.cca() diagnoses linear
# dependencies among constraints/contrasts; VIF > 10 is treated as a warning
# threshold, not an automatic deletion rule.

r2_info <- vegan::RsquareAdj(fit)
vif_values <- vegan::vif.cca(fit)

vif_table <- tibble(
  constraint_or_contrast = names(vif_values),
  vif = as.numeric(vif_values),
  vif_over_10 = as.numeric(vif_values) > 10
)

alias_names <- tryCatch(
  alias(fit, names.only = TRUE),
  error = function(e) character(0)
)

diagnostic_summary <- tibble(
  n_frames = nrow(dat),
  n_dives = n_distinct(dat$dive_id),
  predictors = paste(usable_terms, collapse = " + "),
  r_squared = as.numeric(r2_info$r.squared),
  adjusted_r_squared = as.numeric(r2_info$adj.r.squared),
  max_vif = if (length(vif_values)) max(vif_values, na.rm = TRUE) else NA_real_,
  n_vif_over_10 = sum(vif_values > 10, na.rm = TRUE),
  n_aliased_constraints = length(alias_names),
  permutation_restriction = PERMUTATION_STRATA_COLUMN
)

write_csv(
  diagnostic_summary,
  file.path(out_dir, "24_dbrda_model_diagnostics.csv")
)
write_csv(
  vif_table,
  file.path(out_dir, "24_dbrda_vif.csv")
)
write_lines(
  alias_names,
  file.path(out_dir, "24_dbrda_aliased_constraints.txt")
)

if (any(vif_values > 10, na.rm = TRUE)) {
  warning(
    "dbRDA contains one or more constraint/contrast VIF values > 10. ",
    "Inspect 24_dbrda_vif.csv before ecological interpretation.",
    call. = FALSE
  )
}

perm_control <- permute::how(
  blocks = factor(dat[[PERMUTATION_STRATA_COLUMN]]),
  nperm = N_PERMUTATIONS
)

set.seed(RANDOM_SEED)
overall_test <- anova(
  fit,
  permutations = perm_control
)

set.seed(RANDOM_SEED)
term_test <- anova(
  fit,
  by = "margin",
  permutations = perm_control
)

set.seed(RANDOM_SEED)
axis_test <- anova(
  fit,
  by = "axis",
  permutations = perm_control
)

write_csv(
  add_adjusted_p(overall_test),
  file.path(out_dir, "24_dbrda_overall_test.csv")
)
write_csv(
  add_adjusted_p(term_test),
  file.path(out_dir, "24_dbrda_term_tests.csv")
)
write_csv(
  add_adjusted_p(axis_test),
  file.path(out_dir, "24_dbrda_axis_tests.csv")
)

score_matrix <- scores(fit, display = "sites")
if (is.null(dim(score_matrix)) || ncol(score_matrix) < 2) {
  stop(
    "dbRDA returned fewer than two site-score axes; a two-dimensional CAP figure cannot be drawn.",
    call. = FALSE
  )
}

site_scores <- as.data.frame(score_matrix[, 1:2, drop = FALSE])
site_scores$filename <- rownames(site_scores)
site_scores <- as_tibble(site_scores) %>% left_join(dat, by = "filename")

write_csv(site_scores, file.path(out_dir, "24_dbrda_site_scores.csv"))
saveRDS(fit, file.path(out_dir, "24_dbrda_model.rds"))

axis_names <- names(site_scores)[1:2]

p9 <- ggplot(
  site_scores,
  aes(
    x = .data[[axis_names[1]]],
    y = .data[[axis_names[2]]],
    colour = depth_m
  )
) +
  geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
  scale_colour_viridis_c() +
  coord_equal() +
  labs(
    title = "Constrained community ordination",
    subtitle = paste0(
      "Distance-based RDA using ", DISTANCE_METHOD,
      "; environmental constraints: ", paste(usable_terms, collapse = ", ")
    ),
    x = axis_names[1],
    y = axis_names[2],
    colour = "Depth (m)"
  ) +
  theme_bw(base_size = 11)

save_figure(p9, "Fig_09_CAP_dbRDA_Environment")

write_csv(
  tibble(
    analysis = "CAP/dbRDA environmental model",
    point_sampled_frames = nrow(meta_all),
    total_eligible_community_frames = nrow(meta),
    complete_frames_before_replication_filter =
      n_complete_before_replication_filter,
    complete_frames_used = nrow(dat),
    excluded_for_missing_model_data =
      nrow(meta) - n_complete_before_replication_filter,
    excluded_for_rare_factor_levels =
      n_complete_before_replication_filter - nrow(dat),
    variables = paste(usable_terms, collapse = " + "),
    permutation_design = "blocked_within_dive",
    permutation_block = PERMUTATION_STRATA_COLUMN
  ),
  file.path(out_dir, "24_dbrda_sample_size.csv")
)

# -------------------------------------------------------------------------
# Dive-level constrained ordination for latitude/longitude
# -------------------------------------------------------------------------
# Aggregate ONLY the QC-eligible frame universe, so dive-level and frame-level
# multivariate analyses are based on the same biological inclusion rules.
counts_dive <- a$community_counts %>%
  filter(filename %in% eligible_frames) %>%
  group_by(dive_id, community_unit) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  rename(filename = dive_id)

dive_mat <- community_matrix_from_counts(counts_dive, row_id = "filename")
dive_mat <- dive_mat[, colSums(dive_mat) > 0, drop = FALSE]
dive_mat_t <- force_numeric_matrix(
  transform_community_matrix(dive_mat, COMMUNITY_TRANSFORM),
  "Dive-level dbRDA transformed community matrix"
)

dive_meta <- meta %>%
  group_by(dive_id) %>%
  summarise(
    lat = median_or_na(lat),
    long = median_or_na(long),
    eligible_frames = n(),
    .groups = "drop"
  ) %>%
  filter(
    dive_id %in% rownames(dive_mat_t),
    !is.na(lat),
    !is.na(long)
  )

if (nrow(dive_meta) >= 6 &&
    n_distinct(dive_meta$lat) >= 3 &&
    n_distinct(dive_meta$long) >= 3) {

  dive_comm <- dive_mat_t[dive_meta$dive_id, , drop = FALSE]

  dive_fit <- dbrda(
    dive_comm ~ lat + long,
    data = dive_meta,
    distance = DISTANCE_METHOD
  )

  set.seed(RANDOM_SEED)

    dive_r2_info <- vegan::RsquareAdj(dive_fit)
    dive_vif_values <- vegan::vif.cca(dive_fit)

    write_csv(
      tibble(
        constraint_or_contrast = names(dive_vif_values),
        vif = as.numeric(dive_vif_values),
        vif_over_10 = as.numeric(dive_vif_values) > 10
      ),
      file.path(out_dir, "24_dive_dbrda_latlong_vif.csv")
    )

    write_csv(
      tibble(
        n_dives = nrow(dive_dat),
        predictors = "lat + long",
        r_squared = as.numeric(dive_r2_info$r.squared),
        adjusted_r_squared = as.numeric(dive_r2_info$adj.r.squared),
        max_vif = if (length(dive_vif_values)) max(dive_vif_values, na.rm = TRUE) else NA_real_,
        n_vif_over_10 = sum(dive_vif_values > 10, na.rm = TRUE),
        permutation_restriction = "none_dive_is_unit"
      ),
      file.path(out_dir, "24_dive_dbrda_latlong_diagnostics.csv")
    )

  dive_overall <- anova(
    dive_fit,
    permutations = N_PERMUTATIONS
  )

  set.seed(RANDOM_SEED)
  dive_terms <- anova(
    dive_fit,
    by = "margin",
    permutations = N_PERMUTATIONS
  )

  write_csv(
    add_adjusted_p(dive_overall),
    file.path(out_dir, "24_dive_dbrda_latlong_overall_test.csv")
  )

  write_csv(
    add_adjusted_p(dive_terms),
    file.path(out_dir, "24_dive_dbrda_latlong_term_tests.csv")
  )

  dive_scores_matrix <- scores(dive_fit, display = "sites")

  if (!is.null(dim(dive_scores_matrix)) && ncol(dive_scores_matrix) >= 2) {
    dive_scores <- as.data.frame(
      dive_scores_matrix[, 1:2, drop = FALSE]
    )
    dive_scores$dive_id <- rownames(dive_scores)
    dive_scores <- as_tibble(dive_scores) %>%
      left_join(dive_meta, by = "dive_id")

    write_csv(
      dive_scores,
      file.path(out_dir, "24_dive_dbrda_latlong_site_scores.csv")
    )

    axes <- names(dive_scores)[1:2]

    p10 <- ggplot(
      dive_scores,
      aes(
        x = .data[[axes[1]]],
        y = .data[[axes[2]]],
        colour = lat,
        label = dive_id
      )
    ) +
      geom_point(size = POINT_SIZE + 0.5, alpha = POINT_ALPHA) +
      geom_text(
        check_overlap = TRUE,
        nudge_y = 0.03,
        size = 2.6,
        show.legend = FALSE
      ) +
      coord_equal() +
      labs(
        title = "Dive-level constrained community ordination",
        subtitle = "QC-eligible community data aggregated by dive; constrained by latitude and longitude",
        x = axes[1],
        y = axes[2],
        colour = "Latitude"
      ) +
      theme_bw(base_size = 10)

    p10 <- apply_latitude_colour_scale(p10)
    save_figure(
      p10,
      "Fig_10_CAP_dbRDA_Dive_Latitude",
      width = 10,
      height = 8
    )
  }

  saveRDS(
    dive_fit,
    file.path(out_dir, "24_dive_dbrda_latlong_model.rds")
  )

  write_csv(
    tibble(
      analysis = "Dive-level CAP/dbRDA latitude + longitude",
      eligible_dives = nrow(dive_meta),
      eligible_frames_aggregated = sum(dive_meta$eligible_frames),
      variables = "lat + long",
      permutation_design = "unrestricted_dive_is_unit"
    ),
    file.path(out_dir, "24_dive_dbrda_latlong_sample_size.csv")
  )
}

cat("Point-sampled frames: ", nrow(meta_all), "\n", sep = "")
cat("Community-model eligible frames: ", nrow(meta), "\n", sep = "")
cat("dbRDA complete frames used: ", nrow(dat), "\n", sep = "")
cat("Variables: ", paste(usable_terms, collapse = " + "), "\n", sep = "")
cat("Eligible dives represented: ", n_distinct(meta$dive_id), "\n", sep = "")
