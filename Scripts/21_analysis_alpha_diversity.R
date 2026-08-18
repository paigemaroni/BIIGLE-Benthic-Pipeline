#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 21: alpha diversity\n")
cat("============================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)

if (nrow(alpha) == 0) {
  stop("No eligible biological community data were available for alpha diversity.", call. = FALSE)
}

out_root <- file.path(ANALYSES_DIR, "02_Alpha_Diversity")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_csv(alpha, file.path(out_root, "21_frame_alpha_diversity.csv"))

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    observed_common_id_richness
  ),
  file.path(out_root, "Richness", "21_observed_richness.csv")
)

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    rarefied_common_id_richness
  ),
  file.path(out_root, "Rarefaction", "21_rarefied_richness.csv")
)

write_csv(
  alpha %>% select(filename, dive_id, shannon_common_id_diversity),
  file.path(out_root, "Shannon_Wiener", "21_shannon_wiener.csv")
)

write_csv(
  alpha %>% select(filename, dive_id, inverse_simpson_common_id_diversity),
  file.path(out_root, "Simpson", "21_inverse_simpson.csv")
)

write_csv(
  alpha %>% select(filename, dive_id, pielou_evenness),
  file.path(out_root, "Evenness", "21_pielou_evenness.csv")
)

diversity_summary <- alpha %>%
  summarise(
    frames = n(),
    frames_rarefied = sum(!is.na(rarefied_common_id_richness)),
    median_assigned_biotic_points = median(assigned_biotic_points, na.rm = TRUE),
    median_observed_richness = median(observed_common_id_richness, na.rm = TRUE),
    median_rarefied_richness = median(rarefied_common_id_richness, na.rm = TRUE),
    median_shannon = median(shannon_common_id_diversity, na.rm = TRUE),
    median_inverse_simpson = median(inverse_simpson_common_id_diversity, na.rm = TRUE),
    median_evenness = median(pielou_evenness, na.rm = TRUE)
  )

write_csv(diversity_summary, file.path(out_root, "21_alpha_diversity_summary.csv"))

depth_data <- alpha %>% filter(!is.na(depth_m))

if (nrow(depth_data) > 0) {
  p3 <- ggplot(depth_data, aes(x = depth_m, y = observed_common_id_richness)) +
    geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = "Observed biological-group richness across depth",
      x = "Depth (m; positive downward)",
      y = "Observed common-ID richness"
    ) +
    theme_bw(base_size = 11)
  save_figure(p3, "Fig_03_Richness_Depth")

  p4 <- ggplot(depth_data, aes(x = depth_m, y = shannon_common_id_diversity)) +
    geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = "Shannon diversity across depth",
      x = "Depth (m; positive downward)",
      y = "Shannon-Wiener diversity (H')"
    ) +
    theme_bw(base_size = 11)
  save_figure(p4, "Fig_04_Shannon_Depth")
}

substrate_data <- alpha %>%
  filter(
    !is.na(frame_substrate_class),
    frame_substrate_class != ""
  ) %>%
  add_count(frame_substrate_class, name = "class_n") %>%
  filter(class_n >= MIN_GROUP_N)

if (nrow(substrate_data) > 0 && n_distinct(substrate_data$frame_substrate_class) >= 2) {
  p5 <- ggplot(
    substrate_data,
    aes(
      x = reorder(frame_substrate_class, shannon_common_id_diversity, FUN = median),
      y = shannon_common_id_diversity
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 1.4) +
    coord_flip() +
    labs(
      title = "Shannon diversity by frame substrate class",
      subtitle = paste0("Only substrate classes with at least ", MIN_GROUP_N, " frames shown"),
      x = "Frame substrate class",
      y = "Shannon-Wiener diversity (H')"
    ) +
    theme_bw(base_size = 10)

  save_figure(p5, "Fig_05_Shannon_Substrate", width = 10, height = 9)
}

cat("Alpha-diversity frames: ", format(nrow(alpha), big.mark = ","), "\n", sep = "")
cat("Output: ", file.path(out_root, "21_frame_alpha_diversity.csv"), "\n", sep = "")
