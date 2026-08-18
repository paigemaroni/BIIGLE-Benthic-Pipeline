#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 28: VME occurrence\n")
cat("===========================\n")

ensure_analysis_dirs()
point_a <- prepare_point_analysis()
frames <- point_a$frame_summary
vme <- read_final_vme_data()

check_required_columns(
  vme,
  c(
    "dive_id", "filename", "annotation_id", "top_level", "label_name",
    "common_id_short", "common_id_mid", "common_id_full",
    "depth_m", "temperature", "lat", "long",
    "frame_substrate_class", "frame_relief_class"
  ),
  "vme_annotations.csv"
)

if (any(standardise_top_level(vme$top_level) != "VME", na.rm = TRUE)) {
  stop("vme_annotations.csv contains non-VME rows.", call. = FALSE)
}

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
  left_join(
    observed %>% select(-dive_id),
    by = "filename"
  ) %>%
  mutate(
    has_vme = replace_na(has_vme, FALSE),
    n_vme_annotations = replace_na(n_vme_annotations, 0L),
    n_vme_types = replace_na(n_vme_types, 0L)
  )

out_dir <- file.path(ANALYSES_DIR, "09_VME")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(vme_unique, file.path(out_dir, "28_vme_unique_annotations.csv"))
write_csv(frame_vme, file.path(out_dir, "28_vme_presence_across_point_frames.csv"))

vme_frames_not_in_point_data <- observed %>%
  filter(!filename %in% frames$filename)

write_csv(
  vme_frames_not_in_point_data,
  file.path(out_dir, "28_vme_frames_not_in_point_data.csv")
)

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

p10 <- ggplot(by_dive, aes(x = reorder(dive_id, vme_frame_prevalence), y = vme_frame_prevalence)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "VME frame occurrence by dive",
    x = "Dive",
    y = "Proportion of point-sampled frames with a recorded VME"
  ) +
  theme_bw(base_size = 10)

save_figure(p10, "Fig_11_VME_Occurrence_By_Dive", width = 9, height = 8)

# Exploratory single-predictor logistic models.
predictors <- intersect(
  c(
    "depth_m", "temperature", "lat", "long",
    "frame_substrate_class", "frame_relief_class"
  ),
  names(frame_vme)
)

glm_results <- list()

for (predictor in predictors) {
  dat <- frame_vme %>%
    filter(
      !is.na(.data[[predictor]]),
      as.character(.data[[predictor]]) != ""
    )

  if (nrow(dat) < 15 || n_distinct(dat[[predictor]]) < 2) next

  if (predictor %in% CONTINUOUS_ENVIRONMENTAL_VARIABLES) {
    dat$predictor_value <- safe_num(dat[[predictor]])
    dat <- dat %>% filter(is.finite(predictor_value))
    fit <- tryCatch(
      glm(has_vme ~ predictor_value, data = dat, family = binomial()),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      sm <- summary(fit)$coefficients
      glm_results[[length(glm_results) + 1]] <- tibble(
        predictor = predictor,
        level = NA_character_,
        n = nrow(dat),
        estimate = sm[2, "Estimate"],
        std_error = sm[2, "Std. Error"],
        p_value = sm[2, "Pr(>|z|)"]
      )
    }
  } else {
    dat$predictor_factor <- droplevels(factor(dat[[predictor]]))
    if (nlevels(dat$predictor_factor) < 2) next
    fit <- tryCatch(
      glm(has_vme ~ predictor_factor, data = dat, family = binomial()),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      sm <- summary(fit)$coefficients
      if (nrow(sm) > 1) {
        for (i in 2:nrow(sm)) {
          glm_results[[length(glm_results) + 1]] <- tibble(
            predictor = predictor,
            level = rownames(sm)[i],
            n = nrow(dat),
            estimate = sm[i, "Estimate"],
            std_error = sm[i, "Std. Error"],
            p_value = sm[i, "Pr(>|z|)"]
          )
        }
      }
    }
  }
}

if (length(glm_results) > 0) {
  glm_table <- bind_rows(glm_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = P_ADJUST_METHOD))
  write_csv(glm_table, file.path(out_dir, "28_vme_logistic_models.csv"))
}

cat("Point-sampled frames: ", nrow(frame_vme), "\n", sep = "")
cat("Frames with recorded VME: ", sum(frame_vme$has_vme), "\n", sep = "")
cat("Outputs written to ", out_dir, "\n", sep = "")
