#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 27: focal taxon responses\n")
cat("==================================\n")

a <- prepare_point_analysis()
counts <- a$community_counts
frames <- a$frame_summary
eligible_frame_ids <- frames %>%
  filter(
    n_points >= MODEL_MIN_TOTAL_POINTS,
    n_points <= MODEL_MAX_TOTAL_POINTS,
    n_biotic_assigned >= MIN_BIOTIC_POINTS
  ) %>%
  pull(filename)


out_dir <- file.path(ANALYSES_DIR, "08_Taxon_Responses")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxon_rank <- counts %>%
  group_by(community_unit) %>%
  summarise(
    total_points = sum(abundance),
    frames_present = n_distinct(filename),
    dives_present = n_distinct(dive_id),
    .groups = "drop"
  ) %>%
  arrange(desc(frames_present), desc(total_points), community_unit)

write_csv(taxon_rank, file.path(out_dir, "27_taxon_prevalence.csv"))

top_taxa <- head(taxon_rank$community_unit, TOP_N_TAXA)

grid <- tidyr::expand_grid(
  filename = eligible_frame_ids,
  community_unit = top_taxa
) %>%
  left_join(
    counts %>% select(filename, community_unit, abundance),
    by = c("filename", "community_unit")
  ) %>%
  mutate(abundance = replace_na(abundance, 0)) %>%
  left_join(
    frames %>% select(
      filename, dive_id, n_biotic_assigned,
      depth_m, temperature, lat, long,
      frame_substrate_class, frame_relief_class
    ),
    by = "filename"
  ) %>%
  mutate(
    relative_abundance = if_else(
      n_biotic_assigned > 0,
      abundance / n_biotic_assigned,
      NA_real_
    ),
    presence = abundance > 0
  )

write_csv(grid, file.path(out_dir, "27_top_taxa_frame_data.csv"))

continuous <- intersect(CONTINUOUS_ENVIRONMENTAL_VARIABLES, names(grid))
results <- list()

for (taxon in top_taxa) {
  for (predictor in continuous) {
    dat <- grid %>%
      filter(
        community_unit == taxon,
        !is.na(.data[[predictor]]),
        !is.na(relative_abundance)
      )

    if (nrow(dat) < 10 || n_distinct(dat[[predictor]]) < 3) next

    sp <- safe_spearman(dat$relative_abundance, dat[[predictor]])

    dat$predictor_value <- safe_num(dat[[predictor]])

    glm_fit <- tryCatch(
      glm(
        presence ~ predictor_value,
        data = dat,
        family = binomial()
      ),
      error = function(e) NULL
    )

    if (!is.null(glm_fit)) {
      sm <- summary(glm_fit)$coefficients
      if (nrow(sm) >= 2) {
        results[[length(results) + 1]] <- tibble(
          community_unit = taxon,
          predictor = predictor,
          n = nrow(dat),
          frames_present = sum(dat$presence),
          spearman_rho_relative_abundance = sp$rho,
          spearman_p = sp$p_value,
          logistic_estimate = sm[2, "Estimate"],
          logistic_std_error = sm[2, "Std. Error"],
          logistic_p = sm[2, "Pr(>|z|)"]
        )
      }
    }
  }
}

if (length(results) > 0) {
  result_table <- bind_rows(results) %>%
    mutate(
      spearman_p_adjusted = p.adjust(spearman_p, method = P_ADJUST_METHOD),
      logistic_p_adjusted = p.adjust(logistic_p, method = P_ADJUST_METHOD)
    )
  write_csv(result_table, file.path(out_dir, "27_top_taxon_environment_associations.csv"))
}

cat("Top taxa analysed: ", length(top_taxa), "\n", sep = "")
cat("Outputs written to ", out_dir, "\n", sep = "")
