#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 23: PERMANOVA and PERMDISP\n")
cat("===================================\n")

a <- prepare_point_analysis()
mat <- a$community_matrix
meta <- a$frame_summary

eligible_frames <- community_model_keep(meta, mat)
mat <- mat[eligible_frames, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]
meta <- meta %>% filter(filename %in% rownames(mat))
meta <- meta[match(rownames(mat), meta$filename), , drop = FALSE]
mat_t <- transform_community_matrix(mat, COMMUNITY_TRANSFORM)

out_dir <- file.path(ANALYSES_DIR, "04_PERMANOVA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_tables <- list()
sample_reports <- list()
dispersion_tables <- list()

run_frame_permanova <- function(variable) {
  dat <- meta
  vals <- dat[[variable]]
  keep_rows <- !is.na(vals) & as.character(vals) != ""
  dat <- dat[keep_rows, , drop = FALSE]
  comm <- mat_t[dat$filename, , drop = FALSE]

  if (nrow(dat) < 8) return(NULL)
  if (length(unique(dat[[variable]])) < 2) return(NULL)

  use_strata <- PERMUTATION_STRATA_COLUMN %in% names(dat) &&
    n_distinct(dat[[PERMUTATION_STRATA_COLUMN]]) >= 2

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

  tab <- extract_anova_table(fit) %>%
    mutate(
      model = paste0("frame_", variable),
      predictor = variable,
      sampling_level = "frame",
      permutation_restriction = if (use_strata) PERMUTATION_STRATA_COLUMN else "none",
      n = nrow(dat),
      n_dives = n_distinct(dat$dive_id)
    )

  sample_reports[[length(sample_reports) + 1]] <<- tibble(
    analysis = paste0("PERMANOVA frame: ", variable),
    total_rows_available = nrow(meta),
    rows_complete_for_model = nrow(dat),
    rows_excluded_missing = nrow(meta) - nrow(dat),
    variables = variable
  )

  # PERMDISP for categorical predictors.
  if (variable %in% CATEGORICAL_ENVIRONMENTAL_VARIABLES) {
    group <- factor(dat[[variable]])
    good_groups <- names(which(table(group) >= MIN_GROUP_N))
    keep_group <- group %in% good_groups
    dat_d <- dat[keep_group, , drop = FALSE]
    comm_d <- mat_t[dat_d$filename, , drop = FALSE]
    group_d <- droplevels(factor(dat_d[[variable]]))

    if (nlevels(group_d) >= 2 && nrow(dat_d) >= 2 * MIN_GROUP_N) {
      d <- vegdist(comm_d, method = DISTANCE_METHOD)
      bd <- betadisper(d, group_d, type = "median", bias.adjust = TRUE)
      set.seed(RANDOM_SEED)
      perm <- permutest(bd, permutations = N_PERMUTATIONS)
      dispersion_tables[[length(dispersion_tables) + 1]] <<-
        extract_anova_table(perm$tab) %>%
        mutate(
          predictor = variable,
          n = nrow(dat_d),
          groups = nlevels(group_d)
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

# Combined frame-level environmental model.
combined_terms <- intersect(PERMANOVA_FRAME_TERMS, names(meta))
if (length(combined_terms) > 0) {
  complete <- complete.cases(meta[, combined_terms, drop = FALSE])
  dat <- meta[complete, , drop = FALSE]

  # Drop categorical levels with no remaining observations automatically.
  for (term in intersect(combined_terms, CATEGORICAL_ENVIRONMENTAL_VARIABLES)) {
    dat[[term]] <- droplevels(factor(dat[[term]]))
  }

  comm <- mat_t[dat$filename, , drop = FALSE]

  if (nrow(dat) >= 15 && n_distinct(dat$dive_id) >= 2) {
    formula <- as.formula(
      paste("comm ~", paste(combined_terms, collapse = " + "))
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
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      result_tables[[length(result_tables) + 1]] <-
        extract_anova_table(fit) %>%
        mutate(
          model = "frame_combined_environment",
          predictor = term,
          sampling_level = "frame",
          permutation_restriction = PERMUTATION_STRATA_COLUMN,
          n = nrow(dat),
          n_dives = n_distinct(dat$dive_id)
        )
    }

    sample_reports[[length(sample_reports) + 1]] <- tibble(
      analysis = "PERMANOVA frame combined environment",
      total_rows_available = nrow(meta),
      rows_complete_for_model = nrow(dat),
      rows_excluded_missing = nrow(meta) - nrow(dat),
      variables = paste(combined_terms, collapse = " + ")
    )
  }
}

# Dive-level tests for latitude and longitude.
# This avoids treating repeated frames as independent units for predictors
# that are typically constant within a dive.
counts_dive <- a$community_counts %>%
  group_by(dive_id, community_unit) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

dive_mat <- community_matrix_from_counts(
  counts_dive %>% rename(filename = dive_id),
  row_id = "filename"
)
rownames(dive_mat) <- rownames(dive_mat)
dive_mat_t <- transform_community_matrix(dive_mat, COMMUNITY_TRANSFORM)

dive_meta <- meta %>%
  group_by(dive_id) %>%
  summarise(
    lat = median_or_na(lat),
    long = median_or_na(long),
    .groups = "drop"
  )

for (variable in PERMANOVA_DIVE_TERMS) {
  if (!variable %in% names(dive_meta)) next
  dat <- dive_meta %>% filter(!is.na(.data[[variable]]))
  dat <- dat %>% filter(dive_id %in% rownames(dive_mat_t))
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
    extract_anova_table(fit) %>%
    mutate(
      model = paste0("dive_", variable),
      predictor = variable,
      sampling_level = "dive",
      permutation_restriction = "none",
      n = nrow(dat),
      n_dives = nrow(dat)
    )
}

if (length(result_tables) > 0) {
  write_csv(
    bind_rows(result_tables),
    file.path(out_dir, "23_permanova_results.csv")
  )
}

if (length(dispersion_tables) > 0) {
  write_csv(
    bind_rows(dispersion_tables),
    file.path(out_dir, "23_permdisp_results.csv")
  )
}

if (length(sample_reports) > 0) {
  write_csv(
    bind_rows(sample_reports),
    file.path(out_dir, "23_permanova_sample_sizes.csv")
  )
}

cat("PERMANOVA outputs written to ", out_dir, "\n", sep = "")
