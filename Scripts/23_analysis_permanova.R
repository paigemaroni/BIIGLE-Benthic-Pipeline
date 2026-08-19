#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 23: PERMANOVA and PERMDISP\n")
cat("===================================\n")

a <- prepare_point_analysis()
mat_all <- a$community_matrix
meta_all <- a$frame_summary

eligible_frames <- community_model_keep(meta_all, mat_all)
mat <- mat_all[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]
meta <- meta_all %>% filter(filename %in% rownames(mat))
meta <- meta[match(rownames(mat), meta$filename), , drop = FALSE]

mat_t <- force_numeric_matrix(
  transform_community_matrix(mat, COMMUNITY_TRANSFORM),
  "PERMANOVA transformed community matrix"
)

out_dir <- file.path(ANALYSES_DIR, "04_PERMANOVA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_tables <- list()
sample_reports <- list()
dispersion_tables <- list()

trim_rare_factor_levels <- function(dat, terms, min_n = MIN_GROUP_N) {
  factor_terms <- intersect(terms, CATEGORICAL_ENVIRONMENTAL_VARIABLES)
  if (length(factor_terms) == 0 || nrow(dat) == 0) return(dat)

  # Iterative trimming is needed because removing a rare level from one factor
  # can make another factor level under-replicated in the remaining data.
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

within_dive_variation_count <- function(dat, variable) {
  if (!"dive_id" %in% names(dat) || !variable %in% names(dat)) return(NA_integer_)
  dat %>%
    filter(!is.na(dive_id), !is.na(.data[[variable]])) %>%
    group_by(dive_id) %>%
    summarise(n_values = n_distinct(.data[[variable]]), .groups = "drop") %>%
    summarise(n = sum(n_values >= 2)) %>%
    pull(n)
}

standardise_permanova_table <- function(tab) {
  out <- extract_anova_table(tab)

  if ("Pr(>F)" %in% names(out)) {
    out <- out %>% mutate(p_value = .data[["Pr(>F)"]])
  } else {
    out$p_value <- NA_real_
  }

  if ("R2" %in% names(out)) {
    out <- out %>% mutate(r_squared = .data[["R2"]])
  } else {
    out$r_squared <- NA_real_
  }

  if ("F" %in% names(out)) {
    out <- out %>% mutate(pseudo_f = .data[["F"]])
  } else {
    out$pseudo_f <- NA_real_
  }

  out
}

run_frame_permanova <- function(variable) {
  dat <- meta

  vals <- dat[[variable]]
  keep_rows <- !is.na(vals) & as.character(vals) != ""
  dat <- dat[keep_rows, , drop = FALSE]
  n_before_replication_filter <- nrow(dat)

  if (variable %in% CATEGORICAL_ENVIRONMENTAL_VARIABLES) {
    dat <- trim_rare_factor_levels(dat, variable, MIN_GROUP_N)
    if (nrow(dat) == 0) return(NULL)
    dat[[variable]] <- droplevels(factor(dat[[variable]]))
  }

  if (nrow(dat) < 8) return(NULL)
  if (n_distinct(dat[[variable]]) < 2) return(NULL)

  comm <- mat_t[dat$filename, , drop = FALSE]

  use_strata <- PERMUTATION_STRATA_COLUMN %in% names(dat) &&
    n_distinct(dat[[PERMUTATION_STRATA_COLUMN]]) >= 2

  n_dives_with_variation <- within_dive_variation_count(dat, variable)

  # When inference is restricted within dive, a predictor with no within-dive
  # variation has no valid permutation contrast at this sampling level.
  if (use_strata && is.finite(n_dives_with_variation) &&
      n_dives_with_variation < 1) {
    warning(
      paste0(
        "Skipping frame PERMANOVA for ", variable,
        ": predictor does not vary within any dive."
      )
    )
    return(NULL)
  }

  formula <- as.formula(paste0("comm ~ ", variable))
  set.seed(RANDOM_SEED)

  fit <- adonis2(
    formula,
    data = dat,
    permutations = N_PERMUTATIONS,
    method = DISTANCE_METHOD,
    by = "margin",
    strata = if (use_strata) dat[[PERMUTATION_STRATA_COLUMN]] else NULL
  )

  tab <- standardise_permanova_table(fit) %>%
    mutate(
      model = paste0("frame_", variable),
      predictor = term,
      sampling_level = "frame",
      permutation_restriction = if (use_strata) PERMUTATION_STRATA_COLUMN else "none",
      n = nrow(dat),
      n_dives = n_distinct(dat$dive_id),
      dives_with_within_predictor_variation = n_dives_with_variation,
      minimum_factor_group_n = if (
        variable %in% CATEGORICAL_ENVIRONMENTAL_VARIABLES
      ) MIN_GROUP_N else NA_integer_
    )

  sample_reports[[length(sample_reports) + 1]] <<- tibble(
    analysis = paste0("PERMANOVA frame: ", variable),
    total_point_sampled_frames = nrow(meta_all),
    total_community_eligible_frames = nrow(meta),
    rows_complete_before_replication_filter = n_before_replication_filter,
    rows_used = nrow(dat),
    rows_excluded_missing = nrow(meta) - n_before_replication_filter,
    rows_excluded_rare_factor_levels =
      n_before_replication_filter - nrow(dat),
    variables = variable
  )

  # PERMDISP for categorical predictors, using the SAME replicated levels and
  # a permutation design blocked by dive.
  if (variable %in% CATEGORICAL_ENVIRONMENTAL_VARIABLES) {
    group_d <- droplevels(factor(dat[[variable]]))

    if (nlevels(group_d) >= 2 && nrow(dat) >= 2 * MIN_GROUP_N) {
      d <- vegdist(comm, method = DISTANCE_METHOD)
      bd <- betadisper(d, group_d, type = "median", bias.adjust = TRUE)

      perm_control <- permute::how(
        blocks = factor(dat$dive_id),
        nperm = N_PERMUTATIONS
      )

      set.seed(RANDOM_SEED)
      perm <- permutest(bd, permutations = perm_control)

      dispersion_tables[[length(dispersion_tables) + 1]] <<-
        standardise_permanova_table(perm$tab) %>%
        mutate(
          model = paste0("PERMDISP_", variable),
          predictor = variable,
          sampling_level = "frame",
          permutation_restriction = "dive_id_blocks",
          n = nrow(dat),
          n_dives = n_distinct(dat$dive_id),
          groups = nlevels(group_d),
          minimum_factor_group_n = MIN_GROUP_N
        )
    }
  }

  tab
}

for (variable in PERMANOVA_FRAME_TERMS) {
  if (variable %in% names(meta)) {
    res <- tryCatch(
      run_frame_permanova(variable),
      error = function(e) {
        warning(paste0("Frame PERMANOVA failed for ", variable, ": ", e$message))
        NULL
      }
    )
    if (!is.null(res)) result_tables[[length(result_tables) + 1]] <- res
  }
}

# -------------------------------------------------------------------------
# Combined frame-level environmental model
# -------------------------------------------------------------------------
combined_terms <- intersect(PERMANOVA_FRAME_TERMS, names(meta))

if (length(combined_terms) > 0) {
  complete <- complete.cases(meta[, combined_terms, drop = FALSE])
  dat <- meta[complete, , drop = FALSE]
  n_complete_before_replication_filter <- nrow(dat)

  dat <- trim_rare_factor_levels(dat, combined_terms, MIN_GROUP_N)

  if (nrow(dat) > 0) {
    for (variable in intersect(
      combined_terms,
      CATEGORICAL_ENVIRONMENTAL_VARIABLES
    )) {
      dat[[variable]] <- droplevels(factor(dat[[variable]]))
    }
  }

  usable_terms <- combined_terms[
    vapply(
      combined_terms,
      function(variable) {
        variable %in% names(dat) &&
          nrow(dat) > 0 &&
          n_distinct(dat[[variable]]) >= 2
      },
      logical(1)
    )
  ]

  if (nrow(dat) >= 15 &&
      length(usable_terms) > 0 &&
      n_distinct(dat$dive_id) >= 2) {

    comm <- mat_t[dat$filename, , drop = FALSE]
    formula <- as.formula(
      paste("comm ~", paste(usable_terms, collapse = " + "))
    )

    set.seed(RANDOM_SEED)
    fit <- tryCatch(
      adonis2(
        formula,
        data = dat,
        permutations = N_PERMUTATIONS,
        method = DISTANCE_METHOD,
        by = "margin",
        strata = dat[[PERMUTATION_STRATA_COLUMN]]
      ),
      error = function(e) {
        warning(paste0("Combined frame PERMANOVA failed: ", e$message))
        NULL
      }
    )

    if (!is.null(fit)) {
      result_tables[[length(result_tables) + 1]] <-
        standardise_permanova_table(fit) %>%
        mutate(
          model = "frame_combined_environment",
          predictor = term,
          sampling_level = "frame",
          permutation_restriction = PERMUTATION_STRATA_COLUMN,
          n = nrow(dat),
          n_dives = n_distinct(dat$dive_id),
          dives_with_within_predictor_variation = NA_integer_,
          minimum_factor_group_n = MIN_GROUP_N
        )
    }

    sample_reports[[length(sample_reports) + 1]] <- tibble(
      analysis = "PERMANOVA frame combined environment",
      total_point_sampled_frames = nrow(meta_all),
      total_community_eligible_frames = nrow(meta),
      rows_complete_before_replication_filter =
        n_complete_before_replication_filter,
      rows_used = nrow(dat),
      rows_excluded_missing =
        nrow(meta) - n_complete_before_replication_filter,
      rows_excluded_rare_factor_levels =
        n_complete_before_replication_filter - nrow(dat),
      variables = paste(usable_terms, collapse = " + ")
    )
  }
}

# -------------------------------------------------------------------------
# Dive-level tests for latitude and longitude
# -------------------------------------------------------------------------
# Only community-model-eligible frames contribute to the dive-level community
# matrix. This keeps frame and dive analyses anchored to the SAME QC universe.
counts_dive <- a$community_counts %>%
  filter(filename %in% eligible_frames) %>%
  group_by(dive_id, community_unit) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  rename(filename = dive_id)

dive_mat <- community_matrix_from_counts(counts_dive, row_id = "filename")
dive_mat <- dive_mat[, colSums(dive_mat) > 0, drop = FALSE]
dive_mat_t <- force_numeric_matrix(
  transform_community_matrix(dive_mat, COMMUNITY_TRANSFORM),
  "Dive-level PERMANOVA transformed community matrix"
)

dive_meta <- meta %>%
  group_by(dive_id) %>%
  summarise(
    lat = median_or_na(lat),
    long = median_or_na(long),
    eligible_frames = n(),
    .groups = "drop"
  )

for (variable in PERMANOVA_DIVE_TERMS) {
  if (!variable %in% names(dive_meta)) next

  dat <- dive_meta %>%
    filter(
      !is.na(.data[[variable]]),
      dive_id %in% rownames(dive_mat_t)
    )

  if (nrow(dat) < 6 || n_distinct(dat[[variable]]) < 3) next

  comm <- dive_mat_t[dat$dive_id, , drop = FALSE]
  formula <- as.formula(paste0("comm ~ ", variable))

  set.seed(RANDOM_SEED)
  fit <- adonis2(
    formula,
    data = dat,
    permutations = N_PERMUTATIONS,
    method = DISTANCE_METHOD,
    by = "margin"
  )

  result_tables[[length(result_tables) + 1]] <-
    standardise_permanova_table(fit) %>%
    mutate(
      model = paste0("dive_", variable),
      predictor = term,
      sampling_level = "dive",
      permutation_restriction = "none_dive_is_unit",
      n = nrow(dat),
      n_dives = nrow(dat),
      dives_with_within_predictor_variation = NA_integer_,
      minimum_factor_group_n = NA_integer_
    )

  sample_reports[[length(sample_reports) + 1]] <- tibble(
    analysis = paste0("PERMANOVA dive: ", variable),
    total_point_sampled_frames = nrow(meta_all),
    total_community_eligible_frames = nrow(meta),
    rows_complete_before_replication_filter = nrow(dat),
    rows_used = nrow(dat),
    rows_excluded_missing = nrow(dive_meta) - nrow(dat),
    rows_excluded_rare_factor_levels = 0L,
    variables = variable
  )
}

# -------------------------------------------------------------------------
# Write outputs with BH/FDR adjustment across reported inferential terms
# -------------------------------------------------------------------------
if (length(result_tables) > 0) {
  permanova_out <- bind_rows(result_tables) %>%
    mutate(
      is_inferential_term = !term %in% c("Residual", "Total"),
      p_adjusted_global = if_else(
        is_inferential_term & is.finite(p_value),
        p.adjust(
          ifelse(is_inferential_term, p_value, NA_real_),
          method = P_ADJUST_METHOD
        ),
        NA_real_
      )
    )

  write_csv(
    permanova_out,
    file.path(out_dir, "23_permanova_results.csv")
  )
}

if (length(dispersion_tables) > 0) {
  permdisp_out <- bind_rows(dispersion_tables) %>%
    mutate(
      is_inferential_term = !term %in% c("Residuals", "Residual", "Total"),
      p_adjusted_global = if_else(
        is_inferential_term & is.finite(p_value),
        p.adjust(
          ifelse(is_inferential_term, p_value, NA_real_),
          method = P_ADJUST_METHOD
        ),
        NA_real_
      )
    )

  write_csv(
    permdisp_out,
    file.path(out_dir, "23_permdisp_results.csv")
  )
}

if (length(sample_reports) > 0) {
  write_csv(
    bind_rows(sample_reports),
    file.path(out_dir, "23_permanova_sample_sizes.csv")
  )
}


# -------------------------------------------------------------------------
# Interpretation flags
# -------------------------------------------------------------------------
# This summary does NOT alter the inferential tests or the global BH/FDR
# correction above. It simply makes the intended interpretation explicit.
# In particular, categorical PERMANOVA results are flagged when the matching
# standalone PERMDISP test indicates heterogeneous multivariate dispersion.

if (exists("permanova_out")) {
  inference_summary <- permanova_out %>%
    filter(is_inferential_term) %>%
    transmute(
      model,
      predictor,
      sampling_level,
      permutation_restriction,
      n,
      n_dives,
      p_value,
      p_adjusted_global,
      raw_p_below_alpha = is.finite(p_value) & p_value < ALPHA,
      global_fdr_below_alpha = is.finite(p_adjusted_global) &
        p_adjusted_global < ALPHA
    )

  if (exists("permdisp_out")) {
    dispersion_summary <- permdisp_out %>%
      filter(is_inferential_term) %>%
      transmute(
        predictor,
        permdisp_p_value = p_value,
        permdisp_p_adjusted_global = p_adjusted_global,
        dispersion_differs_global_fdr =
          is.finite(p_adjusted_global) & p_adjusted_global < ALPHA
      )

    inference_summary <- inference_summary %>%
      left_join(dispersion_summary, by = "predictor")
  } else {
    inference_summary <- inference_summary %>%
      mutate(
        permdisp_p_value = NA_real_,
        permdisp_p_adjusted_global = NA_real_,
        dispersion_differs_global_fdr = NA
      )
  }

  inference_summary <- inference_summary %>%
    mutate(
      interpretation_flag = case_when(
        predictor %in% CATEGORICAL_ENVIRONMENTAL_VARIABLES &
          !is.na(dispersion_differs_global_fdr) &
          dispersion_differs_global_fdr ~
          "categorical_result_cautioned_by_dispersion_difference",
        global_fdr_below_alpha ~
          "supported_after_global_fdr",
        raw_p_below_alpha ~
          "raw_signal_not_supported_after_global_fdr",
        TRUE ~
          "not_supported_at_alpha"
      ),
      multiplicity_note =
        "Primary correction = BH/FDR across all reported inferential PERMANOVA terms"
    )

  write_csv(
    inference_summary,
    file.path(out_dir, "23_permanova_interpretation_summary.csv")
  )
}


cat("Point-sampled frames: ", nrow(meta_all), "\n", sep = "")
cat("Community-model eligible frames: ", nrow(meta), "\n", sep = "")
cat("Eligible dives represented: ", n_distinct(meta$dive_id), "\n", sep = "")
cat("PERMANOVA outputs written to ", out_dir, "\n", sep = "")
