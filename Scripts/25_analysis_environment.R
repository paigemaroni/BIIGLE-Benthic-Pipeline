#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 25: environmental associations\n")
cat("=======================================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)
alpha_model <- alpha %>% filter(alpha_model_eligible)

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

frame_predictors <- intersect(
  c("depth_m", "temperature"),
  names(alpha_model)
)

dive_predictors <- intersect(
  c("lat", "long"),
  names(alpha_model)
)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

blocked_within_dive_spearman <- function(dat, nperm = N_PERMUTATIONS) {
  dat <- dat %>%
    group_by(dive_id) %>%
    mutate(
      response_within = response_value - mean(response_value, na.rm = TRUE),
      predictor_within = predictor_value - mean(predictor_value, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(is.finite(response_within), is.finite(predictor_within))

  if (nrow(dat) < 8 ||
      n_distinct(dat$predictor_within) < 3 ||
      n_distinct(dat$dive_id) < 2) {
    return(tibble(
      statistic = NA_real_,
      p_value = NA_real_,
      permutations = nperm
    ))
  }

  obs <- suppressWarnings(
    cor(dat$response_within, dat$predictor_within,
        method = "spearman", use = "complete.obs")
  )

  if (!is.finite(obs)) {
    return(tibble(
      statistic = NA_real_,
      p_value = NA_real_,
      permutations = nperm
    ))
  }

  set.seed(RANDOM_SEED)

  split_idx <- split(seq_len(nrow(dat)), dat$dive_id)

  perm_stats <- replicate(nperm, {
    xperm <- dat$predictor_within
    for (idx in split_idx) {
      if (length(idx) > 1) {
        xperm[idx] <- sample(xperm[idx], length(idx), replace = FALSE)
      }
    }
    suppressWarnings(
      cor(dat$response_within, xperm,
          method = "spearman", use = "complete.obs")
    )
  })

  perm_stats <- perm_stats[is.finite(perm_stats)]

  p <- if (length(perm_stats) == 0) {
    NA_real_
  } else {
    (1 + sum(abs(perm_stats) >= abs(obs))) / (1 + length(perm_stats))
  }

  tibble(
    statistic = obs,
    p_value = p,
    permutations = length(perm_stats)
  )
}

extract_parametric_predictor <- function(fit, term_name) {
  tab <- summary(fit)$p.table
  if (is.null(tab) || !term_name %in% rownames(tab)) {
    return(tibble(
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  tibble(
    estimate = unname(tab[term_name, "Estimate"]),
    std_error = unname(tab[term_name, "Std. Error"]),
    statistic = unname(tab[term_name, grep("value$", colnames(tab))[1]]),
    p_value = unname(tab[term_name, ncol(tab)])
  )
}

# -------------------------------------------------------------------------
# Frame-level depth/temperature associations
# -------------------------------------------------------------------------
# Frames are repeated within dives, so predictor values are centred within
# dive and models include a dive random intercept. This estimates within-dive
# environmental association rather than treating every frame as independent.
# mgcv::s(dive_id, bs="re") is an iid Gaussian random intercept.

frame_results <- list()
frame_gams <- list()
sample_sizes <- list()

for (response in responses) {
  if (!response %in% names(alpha_model)) next

  for (predictor in frame_predictors) {
    dat <- alpha_model %>%
      select(all_of(c(response, predictor, "dive_id"))) %>%
      mutate(
        response_value = safe_num(.data[[response]]),
        predictor_value = safe_num(.data[[predictor]]),
        dive_id = factor(dive_id)
      ) %>%
      filter(
        is.finite(response_value),
        is.finite(predictor_value),
        !is.na(dive_id)
      ) %>%
      group_by(dive_id) %>%
      mutate(
        predictor_dive_mean = mean(predictor_value, na.rm = TRUE),
        predictor_within = predictor_value - predictor_dive_mean
      ) %>%
      ungroup()

    n_dives <- n_distinct(dat$dive_id)
    n_varying_dives <- dat %>%
      group_by(dive_id) %>%
      summarise(n_values = n_distinct(predictor_value), .groups = "drop") %>%
      summarise(n = sum(n_values >= 2)) %>%
      pull(n)

    sample_sizes[[length(sample_sizes) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      sampling_level = "frame",
      effect_scope = "within_dive",
      n = nrow(dat),
      n_dives = n_dives,
      dives_with_within_predictor_variation = n_varying_dives
    )

    if (nrow(dat) < 8 ||
        n_dives < 2 ||
        n_distinct(dat$predictor_within) < 3 ||
        n_varying_dives < 2) next

    # Blocked within-dive Spearman test.
    sp <- blocked_within_dive_spearman(dat)

    frame_results[[length(frame_results) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      analysis_method = "spearman_blocked_within_dive",
      sampling_level = "frame",
      effect_scope = "within_dive",
      n = nrow(dat),
      n_dives = n_dives,
      estimate = sp$statistic,
      std_error = NA_real_,
      statistic = sp$statistic,
      p_value = sp$p_value,
      r_squared = NA_real_,
      deviance_explained = NA_real_,
      edf = NA_real_,
      permutations = sp$permutations
    )

    # Linear predictor + random intercept for dive.
    lin_fit <- tryCatch(
      gam(
        response_value ~ predictor_within + s(dive_id, bs = "re"),
        data = dat,
        method = "REML"
      ),
      error = function(e) NULL
    )

    if (!is.null(lin_fit)) {
      coef_lin <- extract_parametric_predictor(lin_fit, "predictor_within")
      sm <- summary(lin_fit)

      frame_results[[length(frame_results) + 1]] <- tibble(
        response = response,
        predictor = predictor,
        analysis_method = "linear_with_dive_random_intercept",
        sampling_level = "frame",
        effect_scope = "within_dive",
        n = nrow(dat),
        n_dives = n_dives,
        estimate = coef_lin$estimate,
        std_error = coef_lin$std_error,
        statistic = coef_lin$statistic,
        p_value = coef_lin$p_value,
        r_squared = sm$r.sq,
        deviance_explained = sm$dev.expl,
        edf = NA_real_,
        permutations = NA_integer_
      )
    }

    # Log1p response version for non-negative responses.
    if (all(dat$response_value >= 0, na.rm = TRUE)) {
      log_fit <- tryCatch(
        gam(
          log1p(response_value) ~ predictor_within + s(dive_id, bs = "re"),
          data = dat,
          method = "REML"
        ),
        error = function(e) NULL
      )

      if (!is.null(log_fit)) {
        coef_log <- extract_parametric_predictor(log_fit, "predictor_within")
        sm_log <- summary(log_fit)

        frame_results[[length(frame_results) + 1]] <- tibble(
          response = response,
          predictor = predictor,
          analysis_method = "log1p_linear_with_dive_random_intercept",
          sampling_level = "frame",
          effect_scope = "within_dive",
          n = nrow(dat),
          n_dives = n_dives,
          estimate = coef_log$estimate,
          std_error = coef_log$std_error,
          statistic = coef_log$statistic,
          p_value = coef_log$p_value,
          r_squared = sm_log$r.sq,
          deviance_explained = sm_log$dev.expl,
          edf = NA_real_,
          permutations = NA_integer_
        )
      }
    }

    # Smooth within-dive environmental response + random intercept for dive.
    if (nrow(dat) >= 15 &&
        n_distinct(dat$predictor_within) >= 6 &&
        n_varying_dives >= 2) {
      k_value <- min(
        5,
        max(3, n_distinct(dat$predictor_within) - 1)
      )

      form <- as.formula(
        paste0(
          "response_value ~ s(predictor_within, bs='cs', k=",
          k_value,
          ") + s(dive_id, bs='re')"
        )
      )

      fit_gam <- tryCatch(
        gam(form, data = dat, method = "REML"),
        error = function(e) NULL
      )

      if (!is.null(fit_gam)) {
        st <- summary(fit_gam)$s.table
        smooth_row <- grep("^s\\(predictor_within\\)", rownames(st))

        if (length(smooth_row) == 1) {
          sm_gam <- summary(fit_gam)
          frame_gams[[length(frame_gams) + 1]] <- tibble(
            response = response,
            predictor = predictor,
            analysis_method = "gam_with_dive_random_intercept",
            sampling_level = "frame",
            effect_scope = "within_dive",
            n = nrow(dat),
            n_dives = n_dives,
            edf = st[smooth_row, "edf"],
            ref_df = st[smooth_row, "Ref.df"],
            statistic = st[smooth_row, ncol(st) - 1],
            p_value = st[smooth_row, ncol(st)],
            r_squared = sm_gam$r.sq,
            deviance_explained = sm_gam$dev.expl
          )
        }
      }
    }
  }
}

# -------------------------------------------------------------------------
# Dive-level latitude/longitude associations
# -------------------------------------------------------------------------
# Latitude and longitude are effectively constant within dive. Their
# inferential unit is therefore the dive, not the frame.

dive_results <- list()
dive_gams <- list()

for (response in responses) {
  if (!response %in% names(alpha_model)) next

  for (predictor in dive_predictors) {
    dat <- alpha_model %>%
      select(all_of(c(response, predictor, "dive_id"))) %>%
      mutate(
        response_value = safe_num(.data[[response]]),
        predictor_value = safe_num(.data[[predictor]])
      ) %>%
      filter(
        is.finite(response_value),
        is.finite(predictor_value)
      ) %>%
      group_by(dive_id) %>%
      summarise(
        response_value = mean(response_value, na.rm = TRUE),
        predictor_value = median(predictor_value, na.rm = TRUE),
        n_frames = n(),
        .groups = "drop"
      )

    sample_sizes[[length(sample_sizes) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      sampling_level = "dive",
      effect_scope = "between_dive",
      n = nrow(dat),
      n_dives = nrow(dat),
      dives_with_within_predictor_variation = NA_integer_
    )

    if (nrow(dat) < 8 || n_distinct(dat$predictor_value) < 3) next

    sp <- suppressWarnings(
      cor.test(
        dat$response_value,
        dat$predictor_value,
        method = "spearman",
        exact = FALSE
      )
    )

    dive_results[[length(dive_results) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      analysis_method = "spearman_dive_level",
      sampling_level = "dive",
      effect_scope = "between_dive",
      n = nrow(dat),
      n_dives = nrow(dat),
      estimate = unname(sp$estimate),
      std_error = NA_real_,
      statistic = unname(sp$statistic),
      p_value = sp$p.value,
      r_squared = NA_real_,
      deviance_explained = NA_real_,
      edf = NA_real_,
      permutations = NA_integer_
    )

    lm_fit <- lm(response_value ~ predictor_value, data = dat)
    lm_summary <- summary(lm_fit)
    coef_row <- coef(lm_summary)["predictor_value", ]

    dive_results[[length(dive_results) + 1]] <- tibble(
      response = response,
      predictor = predictor,
      analysis_method = "linear_dive_level",
      sampling_level = "dive",
      effect_scope = "between_dive",
      n = nrow(dat),
      n_dives = nrow(dat),
      estimate = coef_row[["Estimate"]],
      std_error = coef_row[["Std. Error"]],
      statistic = coef_row[["t value"]],
      p_value = coef_row[["Pr(>|t|)"]],
      r_squared = lm_summary$r.squared,
      deviance_explained = NA_real_,
      edf = NA_real_,
      permutations = NA_integer_
    )

    if (all(dat$response_value >= 0, na.rm = TRUE)) {
      log_fit <- lm(log1p(response_value) ~ predictor_value, data = dat)
      log_summary <- summary(log_fit)
      coef_log <- coef(log_summary)["predictor_value", ]

      dive_results[[length(dive_results) + 1]] <- tibble(
        response = response,
        predictor = predictor,
        analysis_method = "log1p_linear_dive_level",
        sampling_level = "dive",
        effect_scope = "between_dive",
        n = nrow(dat),
        n_dives = nrow(dat),
        estimate = coef_log[["Estimate"]],
        std_error = coef_log[["Std. Error"]],
        statistic = coef_log[["t value"]],
        p_value = coef_log[["Pr(>|t|)"]],
        r_squared = log_summary$r.squared,
        deviance_explained = NA_real_,
        edf = NA_real_,
        permutations = NA_integer_
      )
    }

    if (nrow(dat) >= 15 && n_distinct(dat$predictor_value) >= 6) {
      k_value <- min(
        5,
        max(3, n_distinct(dat$predictor_value) - 1)
      )
      form <- as.formula(
        paste0(
          "response_value ~ s(predictor_value, bs='cs', k=",
          k_value,
          ")"
        )
      )
      fit_gam <- tryCatch(
        gam(form, data = dat, method = "REML"),
        error = function(e) NULL
      )

      if (!is.null(fit_gam)) {
        st <- summary(fit_gam)$s.table
        sm_gam <- summary(fit_gam)

        dive_gams[[length(dive_gams) + 1]] <- tibble(
          response = response,
          predictor = predictor,
          analysis_method = "gam_dive_level",
          sampling_level = "dive",
          effect_scope = "between_dive",
          n = nrow(dat),
          n_dives = nrow(dat),
          edf = st[1, "edf"],
          ref_df = st[1, "Ref.df"],
          statistic = st[1, ncol(st) - 1],
          p_value = st[1, ncol(st)],
          r_squared = sm_gam$r.sq,
          deviance_explained = sm_gam$dev.expl
        )
      }
    }
  }
}

# -------------------------------------------------------------------------
# Outputs and multiplicity control
# -------------------------------------------------------------------------

regression_tables <- c(frame_results, dive_results)
if (length(regression_tables) > 0) {
  result_table <- bind_rows(regression_tables) %>%
    group_by(analysis_method) %>%
    mutate(
      p_adjusted_within_method = p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_",
        analysis_method
      )
    ) %>%
    ungroup()

  write_csv(
    result_table,
    file.path(out_dir, "25_environmental_regressions.csv")
  )
}

gam_tables <- c(frame_gams, dive_gams)
if (length(gam_tables) > 0) {
  gam_table <- bind_rows(gam_tables) %>%
    group_by(analysis_method) %>%
    mutate(
      p_adjusted_within_method = p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_",
        analysis_method
      )
    ) %>%
    ungroup()

  write_csv(
    gam_table,
    file.path(out_dir, "25_environmental_gams.csv")
  )
}

if (length(sample_sizes) > 0) {
  write_csv(
    bind_rows(sample_sizes),
    file.path(out_dir, "25_environmental_sample_sizes.csv")
  )
}

write_lines(
  c(
    "Analysis 25 hierarchy note",
    "==========================",
    "",
    "depth_m and temperature:",
    "  sampling unit = frame",
    "  inferential scope = within-dive association",
    "  predictor centred within dive",
    "  linear/log-linear/GAM models include dive random intercept",
    "  Spearman significance uses predictor permutations within dive",
    "",
    "lat and long:",
    "  sampling unit = dive",
    "  response aggregated to mean across eligible frames per dive",
    "  coordinate aggregated to median per dive",
    "",
    "Multiplicity:",
    paste0(
      "  ", P_ADJUST_METHOD,
      " correction is applied separately within each analysis method."
    )
  ),
  file.path(out_dir, "25_environmental_analysis_notes.txt")
)

cat("Environmental model outputs written to ", out_dir, "\n", sep = "")
