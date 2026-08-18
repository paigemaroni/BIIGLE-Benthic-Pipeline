#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 24: CAP / distance-based RDA\n")
cat("=====================================\n")

a <- prepare_point_analysis()
mat <- a$community_matrix
meta <- a$frame_summary

eligible_frames <- community_model_keep(meta, mat)
mat <- mat[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]
meta <- meta %>% filter(filename %in% rownames(mat))
meta <- meta[match(rownames(mat), meta$filename), , drop = FALSE]

terms <- intersect(
  c("depth_m", "temperature", "frame_substrate_class", "frame_relief_class"),
  names(meta)
)

complete <- complete.cases(meta[, terms, drop = FALSE])
dat <- meta[complete, , drop = FALSE]
comm <- transform_community_matrix(mat[dat$filename, , drop = FALSE], COMMUNITY_TRANSFORM)

for (term in intersect(terms, CATEGORICAL_ENVIRONMENTAL_VARIABLES)) {
  dat[[term]] <- droplevels(factor(dat[[term]]))
}

if (nrow(dat) < 15 || length(terms) == 0) {
  stop("Too few complete frames for the configured CAP/dbRDA environmental model.", call. = FALSE)
}

out_dir <- file.path(ANALYSES_DIR, "05_CAP_dbRDA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

formula <- as.formula(
  paste("comm ~", paste(terms, collapse = " + "))
)

# vegan currently recommends dbrda() for distance-based RDA.
fit <- dbrda(
  formula,
  data = dat,
  distance = DISTANCE_METHOD
)

set.seed(RANDOM_SEED)
overall_test <- anova(
  fit,
  permutations = N_PERMUTATIONS,
  strata = dat[[PERMUTATION_STRATA_COLUMN]]
)

set.seed(RANDOM_SEED)
term_test <- anova(
  fit,
  by = "margin",
  permutations = N_PERMUTATIONS,
  strata = dat[[PERMUTATION_STRATA_COLUMN]]
)

set.seed(RANDOM_SEED)
axis_test <- anova(
  fit,
  by = "axis",
  permutations = N_PERMUTATIONS,
  strata = dat[[PERMUTATION_STRATA_COLUMN]]
)

write_csv(
  extract_anova_table(overall_test),
  file.path(out_dir, "24_dbrda_overall_test.csv")
)
write_csv(
  extract_anova_table(term_test),
  file.path(out_dir, "24_dbrda_term_tests.csv")
)
write_csv(
  extract_anova_table(axis_test),
  file.path(out_dir, "24_dbrda_axis_tests.csv")
)

score_matrix <- scores(fit, display = "sites")
if (is.null(dim(score_matrix)) || ncol(score_matrix) < 2) {
  stop("dbRDA returned fewer than two site-score axes; a two-dimensional CAP figure cannot be drawn.", call. = FALSE)
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
      "; environmental constraints: ", paste(terms, collapse = ", ")
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
    total_eligible_community_frames = nrow(meta),
    complete_frames_used = nrow(dat),
    excluded_for_missing_model_data = nrow(meta) - nrow(dat),
    variables = paste(terms, collapse = " + "),
    permutation_strata = PERMUTATION_STRATA_COLUMN
  ),
  file.path(out_dir, "24_dbrda_sample_size.csv")
)



# -------------------------------------------------------------------------
# Dive-level constrained ordination for latitude/longitude
# -------------------------------------------------------------------------
counts_dive <- a$community_counts %>%
  group_by(dive_id, community_unit) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  rename(filename = dive_id)

dive_mat <- community_matrix_from_counts(counts_dive, row_id = "filename")
dive_mat <- dive_mat[, colSums(dive_mat) > 0, drop = FALSE]
dive_mat_t <- transform_community_matrix(dive_mat, COMMUNITY_TRANSFORM)

dive_meta <- a$frame_summary %>%
  group_by(dive_id) %>%
  summarise(
    lat = median_or_na(lat),
    long = median_or_na(long),
    .groups = "drop"
  ) %>%
  filter(
    dive_id %in% rownames(dive_mat_t),
    !is.na(lat),
    !is.na(long)
  )

if (nrow(dive_meta) >= 6 && n_distinct(dive_meta$lat) >= 3) {
  dive_comm <- dive_mat_t[dive_meta$dive_id, , drop = FALSE]

  dive_fit <- dbrda(
    dive_comm ~ lat + long,
    data = dive_meta,
    distance = DISTANCE_METHOD
  )

  set.seed(RANDOM_SEED)
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
    extract_anova_table(dive_overall),
    file.path(out_dir, "24_dive_dbrda_latlong_overall_test.csv")
  )

  write_csv(
    extract_anova_table(dive_terms),
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
        subtitle = "Community counts aggregated by dive; constrained by latitude and longitude",
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
}

cat("dbRDA frames used: ", nrow(dat), "\n", sep = "")
cat("Variables: ", paste(terms, collapse = " + "), "\n", sep = "")
