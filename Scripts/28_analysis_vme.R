#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 28: hierarchical VME occurrence\n")
cat("========================================\n")

VME_MIN_PRESENT_FRAMES <- get0("VME_MIN_PRESENT_FRAMES", ifnotfound = 10)
VME_MIN_ABSENT_FRAMES <- get0("VME_MIN_ABSENT_FRAMES", ifnotfound = 10)
VME_MIN_DIVES <- get0("VME_MIN_DIVES", ifnotfound = 3)
VME_MIN_WITHIN_DIVE_RESPONSE_VARIATION <- get0(
  "VME_MIN_WITHIN_DIVE_RESPONSE_VARIATION", ifnotfound = 2
)
VME_MIN_FACTOR_DIVES <- get0("VME_MIN_FACTOR_DIVES", ifnotfound = 2)

ensure_analysis_dirs()
point_a <- prepare_point_analysis()
frames <- point_a$frame_summary
vme <- read_final_vme_data()

check_required_columns(
  vme,
  c(
    "dive_id", "filename", "annotation_id", "top_level", "label_name",
    "common_id_short", "common_id_mid", "common_id_full"
  ),
  "vme_annotations.csv"
)

if (any(standardise_top_level(vme$top_level) != "VME", na.rm = TRUE)) {
  stop("vme_annotations.csv contains non-VME rows.", call. = FALSE)
}

out_dir <- file.path(ANALYSES_DIR, "09_VME")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vme <- vme %>%
  mutate(
    vme_unit = coalesce(
      clean_chr(common_id_full),
      clean_chr(common_id_mid),
      clean_chr(common_id_short),
      clean_chr(label_name)
    ),
    vme_annotation_key = clean_chr(annotation_id)
  )

vme_unique <- vme %>%
  group_by(filename, vme_annotation_key) %>%
  summarise(
    dive_id = first_non_missing(dive_id),
    vme_unit = resolve_single(vme_unit),
    vme_unit_candidates = collapse_unique(vme_unit),
    .groups = "drop"
  )

observed <- vme_unique %>%
  group_by(filename, dive_id) %>%
  summarise(
    has_vme = TRUE,
    n_vme_annotations = n(),
    n_vme_types = n_unique_nonmissing(vme_unit),
    vme_types = collapse_unique(vme_unit),
    .groups = "drop"
  )

frame_vme <- frames %>%
  left_join(observed %>% select(-dive_id), by = "filename") %>%
  mutate(
    has_vme = replace_na(has_vme, FALSE),
    n_vme_annotations = replace_na(n_vme_annotations, 0L),
    n_vme_types = replace_na(n_vme_types, 0L)
  )

write_csv(vme_unique, file.path(out_dir, "28_vme_unique_annotations.csv"))
write_csv(frame_vme, file.path(out_dir, "28_vme_presence_across_point_frames.csv"))

vme_frames_not_in_point_data <- observed %>%
  filter(!filename %in% frames$filename)
write_csv(
  vme_frames_not_in_point_data,
  file.path(out_dir, "28_vme_frames_not_in_point_data.csv")
)

vme_type_prevalence <- vme_unique %>%
  filter(!is.na(vme_unit)) %>%
  group_by(vme_unit) %>%
  summarise(
    unique_annotations = n(),
    frames_present = n_distinct(filename),
    dives_present = n_distinct(dive_id),
    .groups = "drop"
  ) %>%
  arrange(desc(frames_present), desc(unique_annotations))
write_csv(vme_type_prevalence, file.path(out_dir, "28_vme_type_prevalence.csv"))

by_dive <- frame_vme %>%
  group_by(dive_id) %>%
  summarise(
    sampled_frames = n(),
    frames_with_vme = sum(has_vme),
    vme_frame_prevalence = mean(has_vme),
    total_vme_annotations = sum(n_vme_annotations),
    .groups = "drop"
  )
write_csv(by_dive, file.path(out_dir, "28_vme_presence_by_dive.csv"))

p12 <- ggplot(
  by_dive,
  aes(x = reorder(dive_id, vme_frame_prevalence), y = vme_frame_prevalence)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "VME frame occurrence by dive",
    x = "Dive",
    y = "Proportion of point-sampled frames with a recorded VME"
  ) +
  theme_bw(base_size = 10)
save_figure(p12, "Fig_12_VME_Occurrence_By_Dive", width = 9, height = 8)

# -------------------------------------------------------------------------
# Safe GAM helper: warnings and convergence are explicit.
# -------------------------------------------------------------------------
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
        "step failure", "fail(ed|ure)? to converge", "did not converge",
        "iteration limit", "numerically 0 or 1", "not positive definite"
      ),
      collapse = "|"
    )
    critical_warning <- length(warnings) > 0 &&
      any(grepl(critical_pattern, warnings, ignore.case = TRUE))
    fit_converged <- !is.null(fit) && isTRUE(fit$converged)

    list(
      fit = fit,
      warnings = warnings,
      error = error_message,
      critical_warning = critical_warning,
      fit_converged = fit_converged,
      stable = !is.null(fit) && fit_converged && !critical_warning,
      optimizer = paste(optimizer_value, collapse = "/")
    )
  }

  primary <- run_fit(c("outer", "newton"))
  if (is.null(primary$fit) || !primary$fit_converged || primary$critical_warning) {
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

joint_wald_factor <- function(fit, factor_prefix = "group") {
  beta <- coef(fit)
  coef_names <- grep(paste0("^", factor_prefix), names(beta), value = TRUE)
  if (length(coef_names) == 0) {
    return(tibble(df = NA_integer_, statistic = NA_real_, p_value = NA_real_, rank_deficient = TRUE))
  }
  b <- beta[coef_names]
  v <- vcov(fit)[coef_names, coef_names, drop = FALSE]
  keep <- is.finite(b)
  b <- b[keep]
  v <- v[keep, keep, drop = FALSE]
  if (length(b) == 0 || any(!is.finite(v))) {
    return(tibble(df = length(b), statistic = NA_real_, p_value = NA_real_, rank_deficient = TRUE))
  }
  rank_v <- qr(v)$rank
  if (rank_v < length(b)) {
    return(tibble(df = rank_v, statistic = NA_real_, p_value = NA_real_, rank_deficient = TRUE))
  }
  stat <- as.numeric(t(b) %*% solve(v, b))
  tibble(
    df = length(b), statistic = stat,
    p_value = pchisq(stat, df = length(b), lower.tail = FALSE),
    rank_deficient = FALSE
  )
}

fit_diagnostics <- list()
continuous_eligibility <- list()
continuous_results <- list()

# -------------------------------------------------------------------------
# Depth and temperature: frame-level, within-dive effects.
# -------------------------------------------------------------------------
for (predictor in intersect(c("depth_m", "temperature"), names(frame_vme))) {
  dat <- frame_vme %>%
    transmute(
      has_vme,
      dive_id = factor(dive_id),
      predictor_value = safe_num(.data[[predictor]])
    ) %>%
    filter(is.finite(predictor_value), !is.na(dive_id)) %>%
    group_by(dive_id) %>%
    mutate(
      predictor_dive_mean = mean(predictor_value, na.rm = TRUE),
      predictor_within = predictor_value - predictor_dive_mean
    ) %>%
    ungroup()

  n_present <- sum(dat$has_vme)
  n_absent <- sum(!dat$has_vme)
  n_dives <- n_distinct(dat$dive_id)
  response_varying_dives <- dat %>%
    group_by(dive_id) %>%
    summarise(n_response = n_distinct(has_vme), .groups = "drop") %>%
    summarise(n = sum(n_response >= 2)) %>% pull(n)
  predictor_varying_dives <- dat %>%
    group_by(dive_id) %>%
    summarise(n_predictor = n_distinct(predictor_value), .groups = "drop") %>%
    summarise(n = sum(n_predictor >= 3)) %>% pull(n)

  eligible <- n_present >= VME_MIN_PRESENT_FRAMES &&
    n_absent >= VME_MIN_ABSENT_FRAMES &&
    n_dives >= VME_MIN_DIVES &&
    response_varying_dives >= VME_MIN_WITHIN_DIVE_RESPONSE_VARIATION &&
    predictor_varying_dives >= 2

  exclusion_reason <- case_when(
    n_present < VME_MIN_PRESENT_FRAMES ~ "too_few_vme_presence_frames",
    n_absent < VME_MIN_ABSENT_FRAMES ~ "too_few_vme_absence_frames",
    n_dives < VME_MIN_DIVES ~ "too_few_dives",
    response_varying_dives < VME_MIN_WITHIN_DIVE_RESPONSE_VARIATION ~
      "too_few_dives_with_within_dive_vme_variation",
    predictor_varying_dives < 2 ~ "too_few_dives_with_predictor_variation",
    TRUE ~ "eligible"
  )

  continuous_eligibility[[length(continuous_eligibility) + 1]] <- tibble(
    predictor = predictor, n = nrow(dat), n_dives = n_dives,
    frames_with_vme = n_present, frames_without_vme = n_absent,
    dives_with_within_dive_vme_variation = response_varying_dives,
    dives_with_predictor_variation = predictor_varying_dives,
    inferentially_eligible = eligible, exclusion_reason = exclusion_reason
  )

  if (!eligible) next

  fit_info <- fit_gam_safely(
    has_vme ~ predictor_within + s(dive_id, bs = "re"),
    data = dat,
    family = binomial()
  )

  fit_diagnostics[[length(fit_diagnostics) + 1]] <- tibble(
    predictor = predictor,
    response_type = "frame_vme_occurrence",
    analysis_method = "binomial_linear_with_dive_random_intercept",
    optimizer_used = fit_info$optimizer,
    retried_with_bfgs = fit_info$retried_with_bfgs,
    fit_converged = fit_info$fit_converged,
    critical_warning = fit_info$critical_warning,
    stable_for_inference = fit_info$stable,
    warnings = if (length(fit_info$warnings)) paste(unique(fit_info$warnings), collapse = " | ") else NA_character_,
    primary_warnings = if (length(fit_info$primary_warnings)) paste(unique(fit_info$primary_warnings), collapse = " | ") else NA_character_,
    error_message = fit_info$error,
    primary_error = fit_info$primary_error
  )

  if (fit_info$stable) {
    sm <- summary(fit_info$fit)
    tab <- sm$p.table
    continuous_results[[length(continuous_results) + 1]] <- tibble(
      predictor = predictor,
      response_type = "frame_vme_occurrence",
      analysis_method = "binomial_linear_with_dive_random_intercept",
      effect_scope = "within_dive",
      n = nrow(dat), n_dives = n_dives,
      estimate = tab["predictor_within", "Estimate"],
      std_error = tab["predictor_within", "Std. Error"],
      statistic = tab["predictor_within", grep("value$", colnames(tab))[1]],
      p_value = tab["predictor_within", ncol(tab)],
      deviance_explained = sm$dev.expl,
      model_r_squared = sm$r.sq
    )
  }
}

write_csv(bind_rows(continuous_eligibility), file.path(out_dir, "28_continuous_model_eligibility.csv"))
if (length(continuous_results) > 0) {
  continuous_table <- bind_rows(continuous_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(continuous_table, file.path(out_dir, "28_continuous_vme_models.csv"))
}

# -------------------------------------------------------------------------
# Substrate and relief: only within dives where the factor varies.
# -------------------------------------------------------------------------
factor_eligibility <- list()
factor_results <- list()

for (factor_name in intersect(c("frame_substrate_class", "frame_relief_class"), names(frame_vme))) {
  dat0 <- frame_vme %>%
    transmute(
      has_vme,
      dive_id = factor(dive_id),
      group = clean_chr(.data[[factor_name]])
    ) %>%
    filter(!is.na(group), !is.na(dive_id))

  good_levels <- dat0 %>% count(group, name = "n_frames") %>%
    filter(n_frames >= MIN_GROUP_N) %>% pull(group)
  dat1 <- dat0 %>% filter(group %in% good_levels)

  varying_dives <- dat1 %>% group_by(dive_id) %>%
    summarise(n_groups = n_distinct(group), .groups = "drop") %>%
    filter(n_groups >= 2) %>% pull(dive_id)

  dat <- dat1 %>% filter(dive_id %in% varying_dives) %>%
    mutate(dive_id = droplevels(dive_id), group = factor(group))

  final_levels <- dat %>% group_by(group) %>%
    summarise(n_frames = n(), n_dives = n_distinct(dive_id), .groups = "drop") %>%
    filter(n_frames >= MIN_GROUP_N, n_dives >= VME_MIN_FACTOR_DIVES) %>%
    pull(group)
  dat <- dat %>% filter(group %in% final_levels) %>% mutate(group = droplevels(group))

  n_present <- sum(dat$has_vme)
  n_absent <- sum(!dat$has_vme)
  n_groups <- nlevels(dat$group)
  n_dives <- n_distinct(dat$dive_id)

  eligible <- length(varying_dives) >= VME_MIN_FACTOR_DIVES &&
    n_groups >= 2 &&
    n_present >= VME_MIN_PRESENT_FRAMES &&
    n_absent >= VME_MIN_ABSENT_FRAMES

  exclusion_reason <- case_when(
    length(varying_dives) < VME_MIN_FACTOR_DIVES ~ "too_few_dives_with_within_dive_factor_variation",
    n_groups < 2 ~ "fewer_than_two_replicated_factor_levels_after_within_dive_filter",
    n_present < VME_MIN_PRESENT_FRAMES ~ "too_few_vme_presence_frames",
    n_absent < VME_MIN_ABSENT_FRAMES ~ "too_few_vme_absence_frames",
    TRUE ~ "eligible"
  )

  factor_eligibility[[length(factor_eligibility) + 1]] <- tibble(
    predictor = factor_name, n = nrow(dat), n_dives = n_dives,
    factor_levels = n_groups,
    dives_with_within_dive_factor_variation = length(varying_dives),
    frames_with_vme = n_present, frames_without_vme = n_absent,
    inferentially_eligible = eligible, exclusion_reason = exclusion_reason
  )

  if (!eligible) next

  fit_info <- fit_gam_safely(
    has_vme ~ group + s(dive_id, bs = "re"),
    data = dat,
    family = binomial()
  )

  fit_diagnostics[[length(fit_diagnostics) + 1]] <- tibble(
    predictor = factor_name,
    response_type = "frame_vme_occurrence",
    analysis_method = "binomial_factor_with_dive_random_intercept",
    optimizer_used = fit_info$optimizer,
    retried_with_bfgs = fit_info$retried_with_bfgs,
    fit_converged = fit_info$fit_converged,
    critical_warning = fit_info$critical_warning,
    stable_for_inference = fit_info$stable,
    warnings = if (length(fit_info$warnings)) paste(unique(fit_info$warnings), collapse = " | ") else NA_character_,
    primary_warnings = if (length(fit_info$primary_warnings)) paste(unique(fit_info$primary_warnings), collapse = " | ") else NA_character_,
    error_message = fit_info$error,
    primary_error = fit_info$primary_error
  )

  if (fit_info$stable) {
    wald <- joint_wald_factor(fit_info$fit, "group")
    sm <- summary(fit_info$fit)
    factor_results[[length(factor_results) + 1]] <- tibble(
      predictor = factor_name,
      response_type = "frame_vme_occurrence",
      analysis_method = "binomial_factor_with_dive_random_intercept",
      effect_scope = "within_dives_where_factor_varies",
      n = nrow(dat), n_dives = n_dives, factor_levels = n_groups,
      df = wald$df, statistic = wald$statistic, p_value = wald$p_value,
      rank_deficient = wald$rank_deficient,
      deviance_explained = sm$dev.expl,
      model_r_squared = sm$r.sq
    )
  }
}

write_csv(bind_rows(factor_eligibility), file.path(out_dir, "28_factor_model_eligibility.csv"))
if (length(factor_results) > 0) {
  factor_table <- bind_rows(factor_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(factor_table, file.path(out_dir, "28_factor_vme_models.csv"))
}

if (length(fit_diagnostics) > 0) {
  diagnostics <- bind_rows(fit_diagnostics)
  write_csv(diagnostics, file.path(out_dir, "28_model_fit_diagnostics.csv"))
  write_csv(
    diagnostics %>% filter(!stable_for_inference),
    file.path(out_dir, "28_models_excluded_numerical_instability.csv")
  )
}

# -------------------------------------------------------------------------
# Geography: dive is the sampling unit.
# -------------------------------------------------------------------------
dive_geo <- frame_vme %>% group_by(dive_id) %>%
  summarise(
    n_frames = n(),
    vme_prevalence = mean(has_vme),
    lat = median_or_na(lat),
    long = median_or_na(long),
    .groups = "drop"
  )
write_csv(dive_geo, file.path(out_dir, "28_vme_geography_by_dive.csv"))

geo_results <- list()
for (predictor in intersect(c("lat", "long"), names(dive_geo))) {
  dat <- dive_geo %>%
    transmute(
      dive_id,
      vme_prevalence,
      predictor_value = safe_num(.data[[predictor]])
    ) %>%
    filter(is.finite(vme_prevalence), is.finite(predictor_value))

  if (nrow(dat) < 8 || n_distinct(dat$predictor_value) < 3 || n_distinct(dat$vme_prevalence) < 3) next

  sp <- suppressWarnings(cor.test(
    dat$vme_prevalence, dat$predictor_value,
    method = "spearman", exact = FALSE
  ))

  geo_results[[length(geo_results) + 1]] <- tibble(
    predictor = predictor,
    response_type = "dive_vme_frame_prevalence",
    sampling_level = "dive",
    n_dives = nrow(dat),
    rho = unname(sp$estimate),
    p_value = sp$p.value
  )
}

if (length(geo_results) > 0) {
  geo_table <- bind_rows(geo_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(geo_table, file.path(out_dir, "28_geography_vme_associations.csv"))
}

write_lines(
  c(
    "Analysis 28 VME hierarchy and interpretation notes",
    "=================================================",
    "",
    "VME data are kept separate from random-point community data.",
    "A frame with no recorded VME is treated as a VME absence in occurrence models.",
    "This inference is valid only when VME screening/annotation was applied consistently across the analysed frame universe.",
    "",
    "depth_m and temperature:",
    "  sampling unit = frame",
    "  predictor centred within dive",
    "  dive random intercept included",
    "  model warnings/convergence captured explicitly",
    "",
    "substrate and relief:",
    "  inference restricted to dives where the factor varies within dive",
    "  rank-deficient factor tests are flagged rather than interpreted",
    "",
    "latitude and longitude:",
    "  sampling unit = dive",
    "  VME frame prevalence aggregated to one value per dive",
    "",
    paste0("Multiplicity correction: ", P_ADJUST_METHOD, " within each predeclared result family."),
    "",
    "Figure numbering: VME occurrence begins at Fig_12 to avoid overwriting Fig_11 focal-taxon prevalence."
  ),
  file.path(out_dir, "28_analysis_notes.txt")
)

cat("Point-sampled frames: ", nrow(frame_vme), "\n", sep = "")
cat("Frames with recorded VME: ", sum(frame_vme$has_vme), "\n", sep = "")
cat("Unique VME annotations: ", nrow(vme_unique), "\n", sep = "")
cat("VME-only/non-point frames: ", nrow(vme_frames_not_in_point_data), "\n", sep = "")
cat("Outputs written to ", out_dir, "\n", sep = "")
