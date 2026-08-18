#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 26: eligible univariate group comparisons\n")
cat("==================================================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)

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

factors <- intersect(
  c("frame_substrate_class", "frame_relief_class", "dive_id"),
  names(alpha)
)

global_results <- list()
pairwise_results <- list()

for (response in responses) {
  if (!response %in% names(alpha)) next

  for (factor_name in factors) {
    dat <- alpha %>%
      filter(
        !is.na(.data[[response]]),
        !is.na(.data[[factor_name]]),
        as.character(.data[[factor_name]]) != ""
      ) %>%
      mutate(group = factor(.data[[factor_name]])) %>%
      group_by(group) %>%
      filter(n() >= MIN_GROUP_N) %>%
      ungroup() %>%
      mutate(group = droplevels(group))

    n_groups <- nlevels(dat$group)
    if (n_groups < 2) next

    dat$response_value <- safe_num(dat[[response]])

    if (n_groups == 2) {
      lev <- levels(dat$group)
      fit <- t.test(
        dat$response_value[dat$group == lev[[1]]],
        dat$response_value[dat$group == lev[[2]]]
      )
      global_results[[length(global_results) + 1]] <- tibble(
        response = response,
        grouping_variable = factor_name,
        test = "Welch_t_test",
        groups = n_groups,
        n = nrow(dat),
        statistic = unname(fit$statistic),
        df = unname(fit$parameter),
        p_value = fit$p.value
      )
    } else {
      fit <- aov(response_value ~ group, data = dat)
      tab <- summary(fit)[[1]]
      global_results[[length(global_results) + 1]] <- tibble(
        response = response,
        grouping_variable = factor_name,
        test = "ANOVA",
        groups = n_groups,
        n = nrow(dat),
        statistic = tab[["F value"]][1],
        df = tab[["Df"]][1],
        p_value = tab[["Pr(>F)"]][1]
      )

      pw <- pairwise.t.test(
        dat$response_value,
        dat$group,
        p.adjust.method = P_ADJUST_METHOD,
        pool.sd = FALSE
      )

      if (!is.null(pw$p.value)) {
        idx <- which(!is.na(pw$p.value), arr.ind = TRUE)
        if (nrow(idx) > 0) {
          for (i in seq_len(nrow(idx))) {
            r <- rownames(pw$p.value)[idx[i, 1]]
            c <- colnames(pw$p.value)[idx[i, 2]]
            pairwise_results[[length(pairwise_results) + 1]] <- tibble(
              response = response,
              grouping_variable = factor_name,
              group_1 = c,
              group_2 = r,
              p_adjusted = pw$p.value[idx[i, 1], idx[i, 2]]
            )
          }
        }
      }
    }
  }
}

if (length(global_results) > 0) {
  global <- bind_rows(global_results) %>%
    mutate(p_adjusted_across_tests = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(global, file.path(out_dir, "26_global_univariate_tests.csv"))
}

if (length(pairwise_results) > 0) {
  write_csv(
    bind_rows(pairwise_results),
    file.path(out_dir, "26_pairwise_group_comparisons.csv")
  )
}

cat("Only groups with at least ", MIN_GROUP_N, " frames were eligible.\n", sep = "")
cat("Outputs written to ", out_dir, "\n", sep = "")
