#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 20: point-data quality control\n")
cat("=======================================\n")

a <- prepare_point_analysis()
resolved <- a$resolved_points
frames <- a$frame_summary

qc_dir <- file.path(ANALYSES_DIR, "00_Quality_Control")
comp_dir <- file.path(ANALYSES_DIR, "01_Frame_Composition")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(resolved, file.path(qc_dir, "20_resolved_unique_point_annotations.csv"))
write_csv(frames, file.path(qc_dir, "20_frame_qc_summary.csv"))
write_csv(
  frames %>% filter(n_points != TARGET_POINTS_PER_FRAME),
  file.path(qc_dir, "20_frames_not_target_point_count.csv")
)
write_csv(
  resolved %>%
    count(point_class, name = "unique_annotation_points") %>%
    arrange(desc(unique_annotation_points)),
  file.path(qc_dir, "20_point_class_counts.csv")
)

resolution_summary <- resolved %>%
  filter(point_class == "Biotic") %>%
  mutate(
    resolution_status = case_when(
      multiple_community_units ~ "multiple_conflicting_units",
      is.na(community_unit) ~ "unresolved",
      TRUE ~ coalesce(community_unit_source, "resolved_unknown_source")
    )
  ) %>%
  count(resolution_status, name = "biotic_points") %>%
  mutate(percent = 100 * biotic_points / sum(biotic_points))

write_csv(
  resolution_summary,
  file.path(qc_dir, "20_community_unit_resolution_summary.csv")
)

metadata_qc <- tibble(
  variable = c(
    "lat", "long", "depth_m", "temperature",
    "frame_substrate_class", "frame_relief_class"
  ),
  frames_with_value = c(
    sum(!is.na(frames$lat)),
    sum(!is.na(frames$long)),
    sum(!is.na(frames$depth_m)),
    sum(!is.na(frames$temperature)),
    sum(!is.na(frames$frame_substrate_class)),
    sum(!is.na(frames$frame_relief_class))
  )
) %>%
  mutate(
    total_frames = nrow(frames),
    percent_complete = 100 * frames_with_value / total_frames
  )

write_csv(metadata_qc, file.path(qc_dir, "20_frame_metadata_completeness.csv"))

frame_comp <- frames %>%
  group_by(dive_id) %>%
  summarise(
    n_frames = n(),
    mean_pct_biotic = mean(pct_biotic_all, na.rm = TRUE),
    mean_pct_abiotic = mean(pct_abiotic_all, na.rm = TRUE),
    mean_pct_exclude = mean(pct_exclude_all, na.rm = TRUE),
    mean_pct_unsure = mean(pct_unsure_all, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(frame_comp, file.path(comp_dir, "20_frame_composition_by_dive.csv"))

p1 <- ggplot(frames, aes(x = n_points)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  geom_vline(xintercept = TARGET_POINTS_PER_FRAME, linetype = 2) +
  labs(
    title = "Unique random-point annotations per frame",
    subtitle = paste0("Target = ", TARGET_POINTS_PER_FRAME, " unique points per frame"),
    x = "Unique point annotations",
    y = "Number of frames"
  ) +
  theme_bw(base_size = 11)

save_figure(p1, "Fig_01_Point_Count_QC")

long_comp <- frames %>%
  select(dive_id, pct_biotic_all, pct_abiotic_all, pct_exclude_all, pct_unsure_all) %>%
  pivot_longer(
    cols = -dive_id,
    names_to = "class",
    values_to = "percent"
  ) %>%
  mutate(
    class = recode(
      class,
      pct_biotic_all = "Biotic",
      pct_abiotic_all = "Abiotic",
      pct_exclude_all = "Exclude",
      pct_unsure_all = "Unsure"
    )
  )

p2 <- long_comp %>%
  group_by(dive_id, class) %>%
  summarise(mean_percent = mean(percent, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = dive_id, y = mean_percent, fill = class)) +
  geom_col() +
  labs(
    title = "Mean point composition by dive",
    x = "Dive",
    y = "Mean percentage of unique point annotations",
    fill = "Point class"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

save_figure(p2, "Fig_02_Frame_Composition_By_Dive", width = 12, height = 7)

cat("Frames resolved: ", format(nrow(frames), big.mark = ","), "\n", sep = "")
cat(
  "Frames at exactly ", TARGET_POINTS_PER_FRAME, " points: ",
  sum(frames$n_points == TARGET_POINTS_PER_FRAME), "/", nrow(frames), "\n", sep = ""
)
cat("QC outputs: ", qc_dir, "\n", sep = "")
