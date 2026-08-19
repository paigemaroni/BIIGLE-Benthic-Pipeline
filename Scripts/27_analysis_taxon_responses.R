#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 27: hierarchical focal-taxon responses\n")
cat("===============================================\n")

# -------------------------------------------------------------------------
# Configurable defaults
# -------------------------------------------------------------------------
# These fallbacks keep the module runnable with existing repositories. If the
# same names are later added to 00_config.R, those configured values win.

FOCAL_TAXON_MIN_PRESENT_FRAMES <- get0(
  "FOCAL_TAXON_MIN_PRESENT_FRAMES", ifnotfound = 10
)
FOCAL_TAXON_MIN_ABSENT_FRAMES <- get0(
  "FOCAL_TAXON_MIN_ABSENT_FRAMES", ifnotfound = 10
)
FOCAL_TAXON_MIN_DIVES <- get0(
  "FOCAL_TAXON_MIN_DIVES", ifnotfound = 3
)
FOCAL_TAXON_MIN_WITHIN_DIVE_RESPONSE_VARIATION <- get0(
  "FOCAL_TAXON_MIN_WITHIN_DIVE_RESPONSE_VARIATION", ifnotfound = 2
)
FOCAL_TAXON_MIN_FACTOR_DIVES <- get0(
  "FOCAL_TAXON_MIN_FACTOR_DIVES", ifnotfound = 2
)

# -------------------------------------------------------------------------
# Prepare the same QC-eligible universe used by the inferential alpha/community
# analyses. Taxa are ranked within this universe, not across weakly sampled
# frames that would later be excluded from models.
# -------------------------------------------------------------------------

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)
frames <- alpha %>% filter(alpha_model_eligible)
counts <- a$community_counts %>%
  filter(filename %in% frames$filename)

out_dir <- file.path(ANALYSES_DIR, "08_Taxon_Responses")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(frames) == 0 || nrow(counts) == 0) {
  stop("No QC-eligible focal-taxon data are available.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Rank taxa and select the configured number of focal community units.
# -------------------------------------------------------------------------

taxon_rank <- counts %>%
  group_by(community_unit) %>%
  summarise(
    total_points = sum(abundance),
    frames_present = n_distinct(filename),
    dives_present = n_distinct(dive_id),
    .groups = "drop"
  ) %>%
  mutate(
    eligible_frames = nrow(frames),
    eligible_dives = n_distinct(frames$dive_id),
    frame_prevalence_percent = 100 * frames_present / eligible_frames
  ) %>%
  arrange(desc(frames_present), desc(total_points), community_unit)

write_csv(taxon_rank, file.path(out_dir, "27_taxon_prevalence.csv"))

top_taxa <- head(taxon_rank$community_unit, TOP_N_TAXA)

focal_selection <- taxon_rank %>%
  mutate(selected_as_focal_taxon = community_unit %in% top_taxa) %>%
  filter(selected_as_focal_taxon)

write_csv(focal_selection, file.path(out_dir, "27_focal_taxa_selected.csv"))

# -------------------------------------------------------------------------
# Complete frame × focal-taxon grid.
# -------------------------------------------------------------------------

grid <- tidyr::expand_grid(
  filename = frames$filename,
  community_unit = top_taxa
) %>%
  left_join(
    counts %>% select(filename, community_unit, abundance),
    by = c("filename", "community_unit")
  ) %>%
  mutate(abundance = replace_na(abundance, 0L)) %>%
  left_join(
    frames %>% select(
      filename, dive_id, n_points, n_valid, n_biotic_assigned,
      depth_m, temperature, lat, long,
      frame_substrate_class, frame_relief_class
    ),
    by = "filename"
  ) %>%
  mutate(
    presence = abundance > 0,
    taxon_cover_fraction = if_else(
      n_valid > 0,
      abundance / n_valid,
      NA_real_
    ),
    taxon_cover_percent = 100 * taxon_cover_fraction
  )

write_csv(grid, file.path(out_dir, "27_focal_taxon_frame_data.csv"))

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

extract_parametric_term <- function(fit, term_name) {
  tab <- summary(fit)$p.table

  if (is.null(tab) || !term_name %in% rownames(tab)) {
    return(tibble(
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  stat_col <- grep("(z|t) value$", colnames(tab), value = TRUE)
  if (length(stat_col) == 0) stat_col <- colnames(tab)[ncol(tab) - 1]

  tibble(
    estimate = unname(tab[term_name, "Estimate"]),
    std_error = unname(tab[term_name, "Std. Error"]),
    statistic = unname(tab[term_name, stat_col[[1]]]),
    p_value = unname(tab[term_name, ncol(tab)])
  )
}

extract_smooth_term <- function(fit) {
  tab <- summary(fit)$s.table

  if (is.null(tab)) {
    return(tibble(
      edf = NA_real_,
      ref_df = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  idx <- grep("predictor_within", rownames(tab), fixed = TRUE)

  if (length(idx) != 1) {
    return(tibble(
      edf = NA_real_,
      ref_df = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  tibble(
    edf = unname(tab[idx, "edf"]),
    ref_df = unname(tab[idx, "Ref.df"]),
    statistic = unname(tab[idx, ncol(tab) - 1]),
    p_value = unname(tab[idx, ncol(tab)])
  )
}

# Fit mgcv::gam while capturing warnings and checking convergence. If the
# default outer/Newton optimizer emits a step-failure/convergence warning, the
# identical model is retried with outer/BFGS. A p-value is only carried into
# inferential result tables when the final fit is warning-free for critical
# convergence issues and fit$converged is TRUE.
fit_gam_safely <- function(formula, data, family, method = "REML") {
  run_fit <- function(optimizer_value) {
    warnings <- character(0)
    error_message <- NA_character_

    fit <- tryCatch(
      withCallingHandlers(
        gam(
          formula = formula,
          data = data,
          family = family,
          method = method,
          optimizer = optimizer_value
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )

    critical_pattern <- paste(
      c(
        "step failure",
        "fail(ed|ure)? to converge",
        "did not converge",
        "iteration limit",
        "numerically 0 or 1",
        "not positive definite"
      ),
      collapse = "|"
    )

    critical_warning <- if (length(warnings) == 0) {
      FALSE
    } else {
      any(grepl(critical_pattern, warnings, ignore.case = TRUE))
    }

    fit_converged <- !is.null(fit) && isTRUE(fit$converged)
    stable <- !is.null(fit) && fit_converged && !critical_warning

    list(
      fit = fit,
      warnings = warnings,
      error = error_message,
      critical_warning = critical_warning,
      fit_converged = fit_converged,
      stable = stable,
      optimizer = paste(optimizer_value, collapse = "/")
    )
  }

  primary <- run_fit(c("outer", "newton"))

  retry_needed <- (
    is.null(primary$fit) ||
    !primary$fit_converged ||
    primary$critical_warning
  )

  if (retry_needed) {
    retry <- run_fit(c("outer", "bfgs"))
    retry$primary_warnings <- primary$warnings
    retry$primary_error <- primary$error
    retry$retried_with_bfgs <- TRUE
    return(retry)
  }

  primary$primary_warnings <- character(0)
  primary$primary_error <- NA_character_
  primary$retried_with_bfgs <- FALSE
  primary
}

fit_diagnostic_row <- function(
    fit_info, community_unit, predictor, response_type, analysis_method) {
  tibble(
    community_unit = community_unit,
    predictor = predictor,
    response_type = response_type,
    analysis_method = analysis_method,
    optimizer_used = fit_info$optimizer,
    retried_with_bfgs = fit_info$retried_with_bfgs,
    fit_converged = fit_info$fit_converged,
    critical_warning = fit_info$critical_warning,
    stable_for_inference = fit_info$stable,
    warnings = if (length(fit_info$warnings)) {
      paste(unique(fit_info$warnings), collapse = " | ")
    } else {
      NA_character_
    },
    primary_warnings = if (length(fit_info$primary_warnings)) {
      paste(unique(fit_info$primary_warnings), collapse = " | ")
    } else {
      NA_character_
    },
    error_message = fit_info$error,
    primary_error = fit_info$primary_error
  )
}

joint_wald_factor <- function(fit, factor_prefix = "group") {
  beta <- coef(fit)
  coef_names <- grep(paste0("^", factor_prefix), names(beta), value = TRUE)

  if (length(coef_names) == 0) {
    return(tibble(
      df = NA_integer_,
      statistic = NA_real_,
      p_value = NA_real_,
      rank_deficient = TRUE
    ))
  }

  b <- beta[coef_names]
  v <- vcov(fit)[coef_names, coef_names, drop = FALSE]
  keep <- is.finite(b)
  b <- b[keep]
  v <- v[keep, keep, drop = FALSE]

  if (length(b) == 0 || any(!is.finite(v))) {
    return(tibble(
      df = length(b),
      statistic = NA_real_,
      p_value = NA_real_,
      rank_deficient = TRUE
    ))
  }

  rank_v <- qr(v)$rank
  rank_deficient <- rank_v < length(b)

  if (rank_deficient) {
    return(tibble(
      df = rank_v,
      statistic = NA_real_,
      p_value = NA_real_,
      rank_deficient = TRUE
    ))
  }

  stat <- as.numeric(t(b) %*% solve(v, b))

  tibble(
    df = length(b),
    statistic = stat,
    p_value = pchisq(stat, df = length(b), lower.tail = FALSE),
    rank_deficient = FALSE
  )
}

# -------------------------------------------------------------------------
# Continuous depth/temperature models.
# -------------------------------------------------------------------------
# The predictor is centred within dive. This prevents between-dive differences
# from masquerading as a frame-level environmental response. A random intercept
# for dive handles repeated frame observations.
#
# Two response formulations are retained because they answer different
# ecological questions:
#   1. occurrence: was the focal taxon detected in the frame?
#   2. point cover: how many valid random points were assigned to the taxon?
#      The latter is fitted as an overdispersed binomial proportion.

continuous_predictors <- intersect(
  c("depth_m", "temperature"),
  names(grid)
)

continuous_eligibility <- list()
continuous_results <- list()
fit_diagnostics <- list()

for (taxon in top_taxa) {
  for (predictor in continuous_predictors) {
    dat <- grid %>%
      filter(
        community_unit == taxon,
        is.finite(safe_num(.data[[predictor]])),
        n_valid > 0
      ) %>%
      mutate(
        predictor_value = safe_num(.data[[predictor]]),
        dive_id = factor(dive_id)
      ) %>%
      group_by(dive_id) %>%
      mutate(
        predictor_dive_mean = mean(predictor_value, na.rm = TRUE),
        predictor_within = predictor_value - predictor_dive_mean
      ) %>%
      ungroup()

    n_present <- sum(dat$presence)
    n_absent <- sum(!dat$presence)
    n_dives <- n_distinct(dat$dive_id)

    response_varying_dives <- dat %>%
      group_by(dive_id) %>%
      summarise(n_response = n_distinct(presence), .groups = "drop") %>%
      summarise(n = sum(n_response >= 2)) %>%
      pull(n)

    predictor_varying_dives <- dat %>%
      group_by(dive_id) %>%
      summarise(n_predictor = n_distinct(predictor_value), .groups = "drop") %>%
      summarise(n = sum(n_predictor >= 3)) %>%
      pull(n)

    eligible <- (
      n_present >= FOCAL_TAXON_MIN_PRESENT_FRAMES &&
      n_absent >= FOCAL_TAXON_MIN_ABSENT_FRAMES &&
      n_dives >= FOCAL_TAXON_MIN_DIVES &&
      response_varying_dives >=
        FOCAL_TAXON_MIN_WITHIN_DIVE_RESPONSE_VARIATION &&
      predictor_varying_dives >= 2 &&
      n_distinct(dat$predictor_within) >= 4
    )

    exclusion_reason <- case_when(
      n_present < FOCAL_TAXON_MIN_PRESENT_FRAMES ~
        "too_few_presence_frames",
      n_absent < FOCAL_TAXON_MIN_ABSENT_FRAMES ~
        "too_few_absence_frames",
      n_dives < FOCAL_TAXON_MIN_DIVES ~
        "too_few_dives",
      response_varying_dives <
        FOCAL_TAXON_MIN_WITHIN_DIVE_RESPONSE_VARIATION ~
        "too_few_dives_with_within_dive_occurrence_variation",
      predictor_varying_dives < 2 ~
        "too_few_dives_with_predictor_variation",
      n_distinct(dat$predictor_within) < 4 ~
        "insufficient_within_dive_predictor_values",
      TRUE ~ "eligible"
    )

    continuous_eligibility[[length(continuous_eligibility) + 1]] <- tibble(
      community_unit = taxon,
      predictor = predictor,
      n = nrow(dat),
      n_dives = n_dives,
      frames_present = n_present,
      frames_absent = n_absent,
      dives_with_within_dive_occurrence_variation = response_varying_dives,
      dives_with_predictor_variation = predictor_varying_dives,
      inferentially_eligible = eligible,
      exclusion_reason = exclusion_reason
    )

    if (!eligible) next

    # Presence: linear within-dive environmental effect.
    presence_linear_info <- fit_gam_safely(
      presence ~ predictor_within + s(dive_id, bs = "re"),
      data = dat,
      family = binomial()
    )
    fit_diagnostics[[length(fit_diagnostics) + 1]] <- fit_diagnostic_row(
      presence_linear_info, taxon, predictor, "frame_occurrence",
      "binomial_linear_with_dive_random_intercept"
    )
    presence_linear <- presence_linear_info$fit

    if (presence_linear_info$stable) {
      term <- extract_parametric_term(presence_linear, "predictor_within")
      sm <- summary(presence_linear)

      continuous_results[[length(continuous_results) + 1]] <- tibble(
        community_unit = taxon,
        predictor = predictor,
        response_type = "frame_occurrence",
        analysis_method = "binomial_linear_with_dive_random_intercept",
        effect_scope = "within_dive",
        n = nrow(dat),
        n_dives = n_dives,
        frames_present = n_present,
        frames_absent = n_absent,
        estimate = term$estimate,
        std_error = term$std_error,
        edf = NA_real_,
        statistic = term$statistic,
        p_value = term$p_value,
        deviance_explained = sm$dev.expl,
        model_r_squared = sm$r.sq
      )
    }

    # Presence: nonlinear within-dive response.
    if (n_distinct(dat$predictor_within) >= 6) {
      k_value <- min(5, max(3, n_distinct(dat$predictor_within) - 1))
      presence_smooth_info <- fit_gam_safely(
        presence ~ s(predictor_within, bs = "cs", k = k_value) +
          s(dive_id, bs = "re"),
        data = dat,
        family = binomial()
      )
      fit_diagnostics[[length(fit_diagnostics) + 1]] <- fit_diagnostic_row(
        presence_smooth_info, taxon, predictor, "frame_occurrence",
        "binomial_gam_with_dive_random_intercept"
      )
      presence_smooth <- presence_smooth_info$fit

      if (presence_smooth_info$stable) {
        smooth_term <- extract_smooth_term(presence_smooth)
        sm <- summary(presence_smooth)

        continuous_results[[length(continuous_results) + 1]] <- tibble(
          community_unit = taxon,
          predictor = predictor,
          response_type = "frame_occurrence",
          analysis_method = "binomial_gam_with_dive_random_intercept",
          effect_scope = "within_dive",
          n = nrow(dat),
          n_dives = n_dives,
          frames_present = n_present,
          frames_absent = n_absent,
          estimate = NA_real_,
          std_error = NA_real_,
          edf = smooth_term$edf,
          statistic = smooth_term$statistic,
          p_value = smooth_term$p_value,
          deviance_explained = sm$dev.expl,
          model_r_squared = sm$r.sq
        )
      }
    }

    # Point-cover response among valid classified random points.
    cover_linear_info <- fit_gam_safely(
      cbind(abundance, n_valid - abundance) ~
        predictor_within + s(dive_id, bs = "re"),
      data = dat,
      family = quasibinomial()
    )
    fit_diagnostics[[length(fit_diagnostics) + 1]] <- fit_diagnostic_row(
      cover_linear_info, taxon, predictor, "valid_point_cover",
      "quasibinomial_linear_with_dive_random_intercept"
    )
    cover_linear <- cover_linear_info$fit

    if (cover_linear_info$stable) {
      term <- extract_parametric_term(cover_linear, "predictor_within")
      sm <- summary(cover_linear)

      continuous_results[[length(continuous_results) + 1]] <- tibble(
        community_unit = taxon,
        predictor = predictor,
        response_type = "valid_point_cover",
        analysis_method = "quasibinomial_linear_with_dive_random_intercept",
        effect_scope = "within_dive",
        n = nrow(dat),
        n_dives = n_dives,
        frames_present = n_present,
        frames_absent = n_absent,
        estimate = term$estimate,
        std_error = term$std_error,
        edf = NA_real_,
        statistic = term$statistic,
        p_value = term$p_value,
        deviance_explained = sm$dev.expl,
        model_r_squared = sm$r.sq
      )
    }
  }
}

if (length(continuous_eligibility) > 0) {
  write_csv(
    bind_rows(continuous_eligibility),
    file.path(out_dir, "27_continuous_model_eligibility.csv")
  )
}

if (length(continuous_results) > 0) {
  continuous_table <- bind_rows(continuous_results) %>%
    group_by(analysis_method, predictor, response_type) %>%
    mutate(
      p_adjusted_within_model_family =
        p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_", analysis_method,
        "_", predictor,
        "_", response_type
      )
    ) %>%
    ungroup()

  write_csv(
    continuous_table,
    file.path(out_dir, "27_continuous_taxon_models.csv")
  )
}

# -------------------------------------------------------------------------
# Categorical substrate/relief models.
# -------------------------------------------------------------------------
# To avoid a factor merely standing in for dive identity, factor models are
# restricted to dives where the factor itself varies. Levels must also meet
# the configured minimum frame count after this restriction.

factor_predictors <- intersect(
  c("frame_substrate_class", "frame_relief_class"),
  names(grid)
)

factor_eligibility <- list()
factor_results <- list()

for (taxon in top_taxa) {
  for (factor_name in factor_predictors) {
    dat0 <- grid %>%
      filter(
        community_unit == taxon,
        !is.na(.data[[factor_name]]),
        n_valid > 0
      ) %>%
      transmute(
        community_unit,
        filename,
        dive_id = factor(dive_id),
        presence,
        abundance,
        n_valid,
        group = clean_chr(.data[[factor_name]])
      ) %>%
      filter(!is.na(group))

    # Initial level replication filter.
    good_levels <- dat0 %>%
      count(group, name = "n_frames") %>%
      filter(n_frames >= MIN_GROUP_N) %>%
      pull(group)

    dat1 <- dat0 %>% filter(group %in% good_levels)

    varying_dives <- dat1 %>%
      group_by(dive_id) %>%
      summarise(n_groups = n_distinct(group), .groups = "drop") %>%
      filter(n_groups >= 2) %>%
      pull(dive_id)

    # Restrict inference to dives where the factor varies.
    dat <- dat1 %>%
      filter(dive_id %in% varying_dives) %>%
      mutate(
        dive_id = droplevels(dive_id),
        group = factor(group)
      )

    # Re-check level replication after within-dive restriction.
    final_levels <- dat %>%
      group_by(group) %>%
      summarise(
        n_frames = n(),
        n_dives = n_distinct(dive_id),
        .groups = "drop"
      ) %>%
      filter(
        n_frames >= MIN_GROUP_N,
        n_dives >= FOCAL_TAXON_MIN_FACTOR_DIVES
      ) %>%
      pull(group)

    dat <- dat %>%
      filter(group %in% final_levels) %>%
      mutate(group = droplevels(group))

    n_present <- sum(dat$presence)
    n_absent <- sum(!dat$presence)
    n_dives <- n_distinct(dat$dive_id)
    n_groups <- nlevels(dat$group)

    eligible <- (
      n_dives >= FOCAL_TAXON_MIN_FACTOR_DIVES &&
      n_groups >= 2 &&
      n_present >= FOCAL_TAXON_MIN_PRESENT_FRAMES &&
      n_absent >= FOCAL_TAXON_MIN_ABSENT_FRAMES
    )

    exclusion_reason <- case_when(
      length(varying_dives) < FOCAL_TAXON_MIN_FACTOR_DIVES ~
        "too_few_dives_with_within_dive_factor_variation",
      n_groups < 2 ~
        "fewer_than_two_replicated_factor_levels_after_within_dive_filter",
      n_present < FOCAL_TAXON_MIN_PRESENT_FRAMES ~
        "too_few_presence_frames",
      n_absent < FOCAL_TAXON_MIN_ABSENT_FRAMES ~
        "too_few_absence_frames",
      TRUE ~ "eligible"
    )

    factor_eligibility[[length(factor_eligibility) + 1]] <- tibble(
      community_unit = taxon,
      predictor = factor_name,
      n = nrow(dat),
      n_dives = n_dives,
      factor_levels = n_groups,
      dives_with_within_dive_factor_variation = length(varying_dives),
      frames_present = n_present,
      frames_absent = n_absent,
      inferentially_eligible = eligible,
      exclusion_reason = exclusion_reason
    )

    if (!eligible) next

    presence_fit_info <- fit_gam_safely(
      presence ~ group + s(dive_id, bs = "re"),
      data = dat,
      family = binomial()
    )
    fit_diagnostics[[length(fit_diagnostics) + 1]] <- fit_diagnostic_row(
      presence_fit_info, taxon, factor_name, "frame_occurrence",
      "binomial_factor_with_dive_random_intercept"
    )
    presence_fit <- presence_fit_info$fit

    if (presence_fit_info$stable) {
      wald <- joint_wald_factor(presence_fit, "group")
      sm <- summary(presence_fit)

      factor_results[[length(factor_results) + 1]] <- tibble(
        community_unit = taxon,
        predictor = factor_name,
        response_type = "frame_occurrence",
        analysis_method = "binomial_factor_with_dive_random_intercept",
        effect_scope = "within_dives_where_factor_varies",
        n = nrow(dat),
        n_dives = n_dives,
        factor_levels = n_groups,
        df = wald$df,
        statistic = wald$statistic,
        p_value = wald$p_value,
        rank_deficient = wald$rank_deficient,
        deviance_explained = sm$dev.expl,
        model_r_squared = sm$r.sq
      )
    }

    cover_fit_info <- fit_gam_safely(
      cbind(abundance, n_valid - abundance) ~
        group + s(dive_id, bs = "re"),
      data = dat,
      family = quasibinomial()
    )
    fit_diagnostics[[length(fit_diagnostics) + 1]] <- fit_diagnostic_row(
      cover_fit_info, taxon, factor_name, "valid_point_cover",
      "quasibinomial_factor_with_dive_random_intercept"
    )
    cover_fit <- cover_fit_info$fit

    if (cover_fit_info$stable) {
      wald <- joint_wald_factor(cover_fit, "group")
      sm <- summary(cover_fit)

      factor_results[[length(factor_results) + 1]] <- tibble(
        community_unit = taxon,
        predictor = factor_name,
        response_type = "valid_point_cover",
        analysis_method = "quasibinomial_factor_with_dive_random_intercept",
        effect_scope = "within_dives_where_factor_varies",
        n = nrow(dat),
        n_dives = n_dives,
        factor_levels = n_groups,
        df = wald$df,
        statistic = wald$statistic,
        p_value = wald$p_value,
        rank_deficient = wald$rank_deficient,
        deviance_explained = sm$dev.expl,
        model_r_squared = sm$r.sq
      )
    }
  }
}

if (length(factor_eligibility) > 0) {
  write_csv(
    bind_rows(factor_eligibility),
    file.path(out_dir, "27_factor_model_eligibility.csv")
  )
}

if (length(factor_results) > 0) {
  factor_table <- bind_rows(factor_results) %>%
    group_by(analysis_method, predictor, response_type) %>%
    mutate(
      p_adjusted_within_model_family =
        p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_", analysis_method,
        "_", predictor,
        "_", response_type
      )
    ) %>%
    ungroup()

  write_csv(
    factor_table,
    file.path(out_dir, "27_factor_taxon_models.csv")
  )
}

# -------------------------------------------------------------------------
# Fit diagnostics
# -------------------------------------------------------------------------
# Critical convergence warnings are never silently ignored. Models that remain
# unstable after the optimizer retry are omitted from inferential result tables
# and retained here for auditability.

if (length(fit_diagnostics) > 0) {
  fit_diagnostics_table <- bind_rows(fit_diagnostics)
  write_csv(
    fit_diagnostics_table,
    file.path(out_dir, "27_model_fit_diagnostics.csv")
  )

  write_csv(
    fit_diagnostics_table %>% filter(!stable_for_inference),
    file.path(out_dir, "27_models_excluded_numerical_instability.csv")
  )
}

# -------------------------------------------------------------------------
# Dive-level geography.
# -------------------------------------------------------------------------
# Latitude/longitude are effectively dive-level covariates. We therefore
# aggregate taxon occurrence to one prevalence estimate per dive before
# testing geography.

dive_taxon <- grid %>%
  group_by(community_unit, dive_id) %>%
  summarise(
    n_frames = n(),
    frames_present = sum(presence),
    dive_prevalence = mean(presence),
    mean_taxon_cover_percent = mean(taxon_cover_percent, na.rm = TRUE),
    lat = median_or_na(lat),
    long = median_or_na(long),
    .groups = "drop"
  )

write_csv(
  dive_taxon,
  file.path(out_dir, "27_dive_taxon_prevalence.csv")
)

geography_results <- list()

for (taxon in top_taxa) {
  for (predictor in intersect(c("lat", "long"), names(dive_taxon))) {
    dat <- dive_taxon %>%
      filter(
        community_unit == taxon,
        is.finite(safe_num(.data[[predictor]])),
        is.finite(dive_prevalence)
      ) %>%
      mutate(predictor_value = safe_num(.data[[predictor]]))

    if (
      nrow(dat) < 8 ||
      n_distinct(dat$predictor_value) < 3 ||
      n_distinct(dat$dive_prevalence) < 3
    ) next

    sp <- suppressWarnings(
      cor.test(
        dat$dive_prevalence,
        dat$predictor_value,
        method = "spearman",
        exact = FALSE
      )
    )

    geography_results[[length(geography_results) + 1]] <- tibble(
      community_unit = taxon,
      predictor = predictor,
      response_type = "dive_frame_prevalence",
      analysis_method = "spearman_dive_level",
      sampling_level = "dive",
      n_dives = nrow(dat),
      rho = unname(sp$estimate),
      p_value = sp$p.value
    )
  }
}

if (length(geography_results) > 0) {
  geography_table <- bind_rows(geography_results) %>%
    group_by(predictor, response_type) %>%
    mutate(
      p_adjusted_within_geographic_family =
        p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_", predictor,
        "_", response_type
      )
    ) %>%
    ungroup()

  write_csv(
    geography_table,
    file.path(out_dir, "27_geography_taxon_associations.csv")
  )
}

# -------------------------------------------------------------------------
# Prevalence figure.
# -------------------------------------------------------------------------

prevalence_plot <- focal_selection %>%
  mutate(
    community_unit = reorder(community_unit, frame_prevalence_percent)
  ) %>%
  ggplot(aes(x = community_unit, y = frame_prevalence_percent)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Focal community unit",
    y = "Eligible frames present (%)",
    title = "Prevalence of focal benthic community units"
  ) +
  theme_bw()

save_figure(prevalence_plot, "Fig_11_Focal_Taxon_Prevalence")

# -------------------------------------------------------------------------
# Analysis notes.
# -------------------------------------------------------------------------

write_lines(
  c(
    "Analysis 27 hierarchy and interpretation notes",
    "==============================================",
    "",
    paste0("Focal taxa selected: top ", length(top_taxa),
           " community units by prevalence among QC-eligible frames."),
    paste0("QC-eligible frames: ", nrow(frames)),
    paste0("Minimum presence frames for continuous/factor models: ",
           FOCAL_TAXON_MIN_PRESENT_FRAMES),
    paste0("Minimum absence frames for continuous/factor models: ",
           FOCAL_TAXON_MIN_ABSENT_FRAMES),
    "",
    "depth_m and temperature:",
    "  sampling unit = frame",
    "  predictor centred within dive",
    "  dive random intercept included",
    "  occurrence fitted with binomial models",
    "  point cover fitted as taxon points / valid classified points using quasibinomial models",
    "  model R-squared/deviance values describe the whole model, including dive effects",
    "  GAM fits are warning-captured and checked for numerical convergence",
    "  Newton step/convergence failures are retried using outer/BFGS with the identical model",
    "  fits still carrying critical warnings are excluded from inferential p-value tables",
    "",
    "substrate and relief:",
    "  inference restricted to dives where the factor varies within dive",
    paste0("  retained factor levels require at least ", MIN_GROUP_N, " frames"),
    "  factor effect tested jointly using Wald tests from hierarchical GAMs",
    "  rank-deficient factor models are flagged rather than interpreted",
    "",
    "latitude and longitude:",
    "  sampling unit = dive",
    "  focal-taxon frame prevalence aggregated to one value per dive",
    "  Spearman association tested across dives",
    "",
    paste0("Multiplicity correction: ", P_ADJUST_METHOD,
           " within each predeclared model/predictor/response family."),
    "",
    "Interpretation:",
    "  focal taxon models are exploratory taxon-specific follow-ups to the community-level analyses.",
    "  They should not be interpreted as independent species-level causal tests."
  ),
  file.path(out_dir, "27_analysis_notes.txt")
)

cat("QC-eligible frames: ", nrow(frames), "\n", sep = "")
cat("Focal taxa selected: ", length(top_taxa), "\n", sep = "")
cat("Outputs written to ", out_dir, "\n", sep = "")
