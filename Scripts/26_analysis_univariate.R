#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 26: hierarchical univariate group comparisons\n")
cat("======================================================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)
alpha_model <- alpha %>% filter(alpha_model_eligible)

out_dir <- file.path(ANALYSES_DIR, "07_Univariate_Tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

responses <- c(
  "observed_common_id_richness",
  "rarefied_common_id_richness",
  "shannon_common_id_diversity",
  "inverse_simpson_common_id_diversity",
  "pielou_evenness",
  "pct_biotic_all"
)

# Substrate and relief vary at frame level but frames are nested within dives.
# dive_id itself is NOT tested using frames as independent replicates.
frame_factors <- intersect(
  c("frame_substrate_class", "frame_relief_class"),
  names(alpha_model)
)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

blocked_partial_f_test <- function(dat, nperm = N_PERMUTATIONS) {
  dat <- dat %>%
    mutate(
      dive_id = factor(dive_id),
      group = droplevels(factor(group))
    )

  if (nrow(dat) < 8 ||
      nlevels(dat$group) < 2 ||
      n_distinct(dat$dive_id) < 2) {
    return(tibble(
      df_factor = NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      partial_r_squared = NA_real_,
      permutations = 0L,
      aliased_coefficients = NA_integer_
    ))
  }

  reduced <- lm(response_value ~ dive_id, data = dat)
  full <- lm(response_value ~ dive_id + group, data = dat)

  cmp <- anova(reduced, full)

  if (nrow(cmp) < 2 || !is.finite(cmp$F[2])) {
    return(tibble(
      df_factor = if (nrow(cmp) >= 2) cmp$Df[2] else NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      partial_r_squared = NA_real_,
      permutations = 0L,
      aliased_coefficients = sum(is.na(coef(full)))
    ))
  }

  observed_f <- unname(cmp$F[2])
  rss_reduced <- unname(cmp$RSS[1])
  rss_full <- unname(cmp$RSS[2])

  partial_r2 <- if (
    is.finite(rss_reduced) &&
    rss_reduced > 0 &&
    is.finite(rss_full)
  ) {
    (rss_reduced - rss_full) / rss_reduced
  } else {
    NA_real_
  }

  split_idx <- split(seq_len(nrow(dat)), dat$dive_id)
  set.seed(RANDOM_SEED)

  perm_f <- replicate(nperm, {
    perm_group <- as.character(dat$group)

    for (idx in split_idx) {
      if (length(idx) > 1 && n_distinct(perm_group[idx]) > 1) {
        perm_group[idx] <- sample(
          perm_group[idx],
          length(idx),
          replace = FALSE
        )
      }
    }

    perm_dat <- dat
    perm_dat$group_perm <- factor(
      perm_group,
      levels = levels(dat$group)
    )

    perm_full <- tryCatch(
      lm(response_value ~ dive_id + group_perm, data = perm_dat),
      error = function(e) NULL
    )

    if (is.null(perm_full)) return(NA_real_)

    perm_cmp <- tryCatch(
      anova(reduced, perm_full),
      error = function(e) NULL
    )

    if (is.null(perm_cmp) ||
        nrow(perm_cmp) < 2 ||
        !is.finite(perm_cmp$F[2])) {
      return(NA_real_)
    }

    unname(perm_cmp$F[2])
  })

  perm_f <- perm_f[is.finite(perm_f)]

  p <- if (length(perm_f) == 0) {
    NA_real_
  } else {
    (1 + sum(perm_f >= observed_f)) / (1 + length(perm_f))
  }

  tibble(
    df_factor = unname(cmp$Df[2]),
    f_statistic = observed_f,
    p_value = p,
    partial_r_squared = partial_r2,
    permutations = length(perm_f),
    aliased_coefficients = sum(is.na(coef(full)))
  )
}

prepare_factor_data <- function(response, factor_name) {
  dat <- alpha_model %>%
    transmute(
      filename,
      dive_id = factor(dive_id),
      response_value = safe_num(.data[[response]]),
      group = clean_chr(.data[[factor_name]])
    ) %>%
    filter(
      is.finite(response_value),
      !is.na(dive_id),
      !is.na(group)
    )

  before_n <- nrow(dat)

  eligible_levels <- dat %>%
    count(group, name = "n_frames") %>%
    filter(n_frames >= MIN_GROUP_N) %>%
    pull(group)

  dat <- dat %>%
    filter(group %in% eligible_levels) %>%
    mutate(group = factor(group))

  variation <- dat %>%
    group_by(dive_id) %>%
    summarise(
      n_groups = n_distinct(group),
      .groups = "drop"
    )

  list(
    data = dat,
    rows_before_level_filter = before_n,
    rows_after_level_filter = nrow(dat),
    n_groups = nlevels(dat$group),
    n_dives = n_distinct(dat$dive_id),
    dives_with_within_factor_variation = sum(variation$n_groups >= 2)
  )
}

# -------------------------------------------------------------------------
# Global within-dive factor tests
# -------------------------------------------------------------------------

global_results <- list()
eligibility_results <- list()
pairwise_results <- list()
pairwise_eligibility <- list()

for (response in responses) {
  if (!response %in% names(alpha_model)) next

  for (factor_name in frame_factors) {
    prep <- prepare_factor_data(response, factor_name)
    dat <- prep$data

    eligible_for_inference <- (
      prep$n_groups >= 2 &&
      prep$n_dives >= 2 &&
      prep$dives_with_within_factor_variation >= 2
    )

    exclusion_reason <- case_when(
      prep$n_groups < 2 ~ "fewer_than_two_replicated_factor_levels",
      prep$n_dives < 2 ~ "fewer_than_two_dives",
      prep$dives_with_within_factor_variation < 2 ~
        "fewer_than_two_dives_with_within_factor_variation",
      TRUE ~ "eligible"
    )

    eligibility_results[[length(eligibility_results) + 1]] <- tibble(
      response = response,
      grouping_variable = factor_name,
      total_alpha_model_frames = nrow(alpha_model),
      rows_complete_before_replication_filter =
        prep$rows_before_level_filter,
      rows_used_after_replication_filter =
        prep$rows_after_level_filter,
      groups_retained = prep$n_groups,
      n_dives = prep$n_dives,
      dives_with_within_factor_variation =
        prep$dives_with_within_factor_variation,
      minimum_group_n = MIN_GROUP_N,
      inferentially_eligible = eligible_for_inference,
      exclusion_reason = exclusion_reason
    )

    if (!eligible_for_inference) next

    test <- blocked_partial_f_test(dat)

    global_results[[length(global_results) + 1]] <- tibble(
      response = response,
      grouping_variable = factor_name,
      test = "blocked_partial_F_permutation",
      effect_scope = "within_dive",
      n = nrow(dat),
      n_dives = n_distinct(dat$dive_id),
      groups = nlevels(dat$group),
      dives_with_within_factor_variation =
        prep$dives_with_within_factor_variation,
      df_factor = test$df_factor,
      statistic = test$f_statistic,
      p_value = test$p_value,
      partial_r_squared = test$partial_r_squared,
      permutations = test$permutations,
      aliased_coefficients = test$aliased_coefficients,
      permutation_restriction = "factor_labels_permuted_within_dive"
    )

    # ---------------------------------------------------------------------
    # Pairwise tests
    # ---------------------------------------------------------------------
    levs <- levels(dat$group)
    pairs <- combn(levs, 2, simplify = FALSE)

    for (pair in pairs) {
      pair_dat <- dat %>%
        filter(as.character(group) %in% pair) %>%
        mutate(
          group = factor(
            as.character(group),
            levels = pair
          )
        )

      cooccurrence <- pair_dat %>%
        group_by(dive_id) %>%
        summarise(
          n_groups = n_distinct(group),
          .groups = "drop"
        ) %>%
        filter(n_groups == 2)

      cooccurring_dives <- as.character(cooccurrence$dive_id)

      pair_dat <- pair_dat %>%
        filter(as.character(dive_id) %in% cooccurring_dives) %>%
        mutate(
          dive_id = droplevels(dive_id),
          group = droplevels(group)
        )

      group_counts <- pair_dat %>%
        count(group, name = "n_frames")

      pair_eligible <- (
        length(cooccurring_dives) >= 2 &&
        nlevels(pair_dat$group) == 2 &&
        nrow(group_counts) == 2 &&
        all(group_counts$n_frames >= MIN_GROUP_N)
      )

      pair_reason <- case_when(
        length(cooccurring_dives) < 2 ~
          "fewer_than_two_dives_where_both_groups_cooccur",
        nlevels(pair_dat$group) < 2 ~
          "one_group_absent_after_cooccurrence_filter",
        nrow(group_counts) < 2 ||
          any(group_counts$n_frames < MIN_GROUP_N) ~
          "minimum_group_frame_count_not_met_after_cooccurrence_filter",
        TRUE ~ "eligible"
      )

      pairwise_eligibility[[length(pairwise_eligibility) + 1]] <- tibble(
        response = response,
        grouping_variable = factor_name,
        group_1 = pair[[1]],
        group_2 = pair[[2]],
        cooccurring_dives = length(cooccurring_dives),
        n = nrow(pair_dat),
        minimum_group_n = MIN_GROUP_N,
        inferentially_eligible = pair_eligible,
        exclusion_reason = pair_reason
      )

      if (!pair_eligible) next

      pair_test <- blocked_partial_f_test(pair_dat)

      full <- lm(response_value ~ dive_id + group, data = pair_dat)
      coef_names <- grep("^group", names(coef(full)), value = TRUE)
      adjusted_difference_group2_minus_group1 <- if (length(coef_names) == 1) {
        unname(coef(full)[coef_names])
      } else {
        NA_real_
      }

      means <- pair_dat %>%
        group_by(group) %>%
        summarise(
          raw_mean = mean(response_value, na.rm = TRUE),
          raw_median = median(response_value, na.rm = TRUE),
          n_frames = n(),
          .groups = "drop"
        )

      mean_1 <- means$raw_mean[means$group == pair[[1]]]
      mean_2 <- means$raw_mean[means$group == pair[[2]]]

      pairwise_results[[length(pairwise_results) + 1]] <- tibble(
        response = response,
        grouping_variable = factor_name,
        group_1 = pair[[1]],
        group_2 = pair[[2]],
        test = "blocked_pairwise_partial_F_permutation",
        effect_scope = "within_dive_cooccurring_groups_only",
        n = nrow(pair_dat),
        n_dives = n_distinct(pair_dat$dive_id),
        raw_mean_group_1 = ifelse(length(mean_1), mean_1, NA_real_),
        raw_mean_group_2 = ifelse(length(mean_2), mean_2, NA_real_),
        raw_mean_difference_group2_minus_group1 =
          ifelse(length(mean_1) && length(mean_2),
                 mean_2 - mean_1, NA_real_),
        dive_adjusted_difference_group2_minus_group1 =
          adjusted_difference_group2_minus_group1,
        statistic = pair_test$f_statistic,
        p_value = pair_test$p_value,
        partial_r_squared = pair_test$partial_r_squared,
        permutations = pair_test$permutations,
        permutation_restriction =
          "factor_labels_permuted_within_cooccurring_dives"
      )
    }
  }
}

# -------------------------------------------------------------------------
# Dive-level descriptive summaries
# -------------------------------------------------------------------------
# A test of response ~ dive_id using individual frames would use frames as
# pseudoreplicates of each dive. Therefore dive_id is summarised descriptively
# here but is not assigned an inferential p-value.

dive_descriptive <- bind_rows(lapply(responses, function(response) {
  if (!response %in% names(alpha_model)) return(NULL)

  alpha_model %>%
    transmute(
      dive_id,
      response_value = safe_num(.data[[response]])
    ) %>%
    filter(is.finite(response_value), !is.na(dive_id)) %>%
    group_by(dive_id) %>%
    summarise(
      response = response,
      n_frames = n(),
      mean = mean(response_value),
      sd = sd(response_value),
      median = median(response_value),
      q25 = quantile(response_value, 0.25),
      q75 = quantile(response_value, 0.75),
      .groups = "drop"
    ) %>%
    select(response, everything())
}))

# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------

if (length(eligibility_results) > 0) {
  write_csv(
    bind_rows(eligibility_results),
    file.path(out_dir, "26_factor_test_eligibility.csv")
  )
}

if (length(global_results) > 0) {
  global <- bind_rows(global_results) %>%
    mutate(
      p_adjusted_across_global_factor_tests =
        p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_across_all_eligible_global_factor_tests"
      )
    )

  write_csv(
    global,
    file.path(out_dir, "26_global_hierarchical_factor_tests.csv")
  )
}

if (length(pairwise_eligibility) > 0) {
  write_csv(
    bind_rows(pairwise_eligibility),
    file.path(out_dir, "26_pairwise_test_eligibility.csv")
  )
}

if (length(pairwise_results) > 0) {
  pairwise <- bind_rows(pairwise_results) %>%
    group_by(response, grouping_variable) %>%
    mutate(
      p_adjusted_within_response_factor =
        p.adjust(p_value, method = P_ADJUST_METHOD),
      p_adjust_family = paste0(
        P_ADJUST_METHOD,
        "_within_response_and_factor"
      )
    ) %>%
    ungroup()

  write_csv(
    pairwise,
    file.path(out_dir, "26_pairwise_hierarchical_factor_tests.csv")
  )
}

write_csv(
  dive_descriptive,
  file.path(out_dir, "26_dive_descriptive_summaries.csv")
)

write_lines(
  c(
    "Analysis 26 hierarchy and test eligibility",
    "==========================================",
    "",
    "Why ordinary Welch t-tests / ANOVA are not used here:",
    "  The candidate frame-level grouping variables (substrate and relief)",
    "  are observed repeatedly within dives. Individual frames therefore",
    "  are not independent replicates of those groups.",
    "",
    "Primary factor test:",
    "  partial F for grouping variable after controlling for dive as a",
    "  fixed blocking factor; significance from factor-label permutations",
    "  restricted within dive.",
    "",
    "Global factor inference requires:",
    paste0("  at least ", MIN_GROUP_N, " frames per retained factor level"),
    "  at least two retained levels",
    "  at least two dives",
    "  at least two dives showing within-dive factor variation",
    "",
    "Pairwise inference requires:",
    "  both groups to co-occur within at least two dives",
    paste0("  at least ", MIN_GROUP_N,
           " frames per group after restricting to co-occurring dives"),
    "",
    "Dive comparisons:",
    "  descriptive only; frames are not treated as independent replicates",
    "  of a dive-level factor.",
    "",
    paste0(
      "Multiplicity: ",
      P_ADJUST_METHOD,
      " is applied across eligible global tests and separately within",
      " each response × factor family for pairwise tests."
    )
  ),
  file.path(out_dir, "26_analysis_notes.txt")
)

cat("Alpha-model frames: ", nrow(alpha_model), "\n", sep = "")
cat("Frame-level factors considered: ",
    paste(frame_factors, collapse = ", "), "\n", sep = "")
cat("Classical Welch/ANOVA tests: not inferentially eligible under nested sampling.\n")
cat("Hierarchical outputs written to ", out_dir, "\n", sep = "")
