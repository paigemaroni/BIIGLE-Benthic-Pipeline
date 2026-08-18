#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 25: environmental associations\n")
cat("=======================================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)

out_dir <- file.path(ANALYSES_DIR, "06_Environmental_Associations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

responses <- c(
  "observed_common_id_richness",
  "rarefied_common_id_richness",
  "shannon_common_id_diversity",
  "inverse_simpson_common_id_diversity",
  "pielou_evenness",
  "pct_biotic_all"
)

predictors <- intersect(
  CONTINUOUS_ENVIRONMENTAL_VARIABLES,
  names(alpha)
)

results <- list()
gam_results <- list()

for (response in responses) {
  if (!response %in% names(alpha)) next

  for (predictor in predictors) {
    dat <- alpha %>%
      select(all_of(c(response, predictor, "dive_id"))) %>%
      filter(
        is.finite(safe_num(.data[[response]])),
        is.finite(safe_num(.data[[predictor]]))
      ) %>%
      mutate(
        response_value = safe_num(.data[[response]]),
        predictor_value = safe_num(.data[[predictor]])
      )

    if (nrow(dat) < 8 || n_distinct(dat$predictor_value) < 3) next

    sp <- safe_spearman(dat$response_value, dat$predictor_value)

    lm_fit <- lm(response_value ~ predictor_value, data = dat)
    lm_summary <- summary(lm_fit)
    coef_row <- coef(lm_summary)["predictor_value", ]

    results[[length(results) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      model = "linear",
      n = nrow(dat),
      estimate = coef_row[["Estimate"]],
      std_error = coef_row[["Std. Error"]],
      statistic = coef_row[["t value"]],
      p_value = coef_row[["Pr(>|t|)"]],
      r_squared = lm_summary$r.squared,
      spearman_rho = sp$rho,
      spearman_p = sp$p_value
    )

    # Log1p response model only for non-negative responses.
    if (all(dat$response_value >= 0, na.rm = TRUE)) {
      log_fit <- lm(log1p(response_value) ~ predictor_value, data = dat)
      log_summary <- summary(log_fit)
      coef_log <- coef(log_summary)["predictor_value", ]

      results[[length(results) + 1]] <- tibble(
        response = response,
        predictor = predictor,
        model = "log1p_linear",
        n = nrow(dat),
        estimate = coef_log[["Estimate"]],
        std_error = coef_log[["Std. Error"]],
        statistic = coef_log[["t value"]],
        p_value = coef_log[["Pr(>|t|)"]],
        r_squared = log_summary$r.squared,
        spearman_rho = sp$rho,
        spearman_p = sp$p_value
      )
    }

    # Exploratory GAM where replication/gradient support it.
    if (nrow(dat) >= 15 && n_distinct(dat$predictor_value) >= 6) {
      k_value <- min(5, max(3, n_distinct(dat$predictor_value) - 1))
      form <- as.formula(
        paste0("response_value ~ s(predictor_value, bs='cs', k=", k_value, ")")
      )
      fit_gam <- tryCatch(
        gam(form, data = dat, method = "REML"),
        error = function(e) NULL
      )
      if (!is.null(fit_gam)) {
        st <- summary(fit_gam)$s.table
        gam_results[[length(gam_results) + 1]] <- tibble(
          response = response,
          predictor = predictor,
          n = nrow(dat),
          edf = st[1, "edf"],
          ref_df = st[1, "Ref.df"],
          statistic = st[1, ncol(st) - 1],
          p_value = st[1, ncol(st)],
          deviance_explained = summary(fit_gam)$dev.expl
        )
      }
    }
  }
}

if (length(results) > 0) {
  result_table <- bind_rows(results) %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD),
      spearman_p_adjusted = p.adjust(spearman_p, method = P_ADJUST_METHOD)
    )
  write_csv(result_table, file.path(out_dir, "25_environmental_regressions.csv"))
}

if (length(gam_results) > 0) {
  gam_table <- bind_rows(gam_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(gam_table, file.path(out_dir, "25_environmental_gams.csv"))
}

cat("Environmental model outputs written to ", out_dir, "\n", sep = "")
