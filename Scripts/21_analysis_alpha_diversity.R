#!/usr/bin/env Rscript
source("Scripts/00_config.R")
source("Scripts/01_analysis_helpers.R")

cat("Analysis 21: alpha diversity\n")
cat("============================\n")

a <- prepare_point_analysis()
alpha <- frame_alpha_diversity(a$community_matrix, a$frame_summary)

if (nrow(alpha) == 0) {
  stop("No point-sampled frame data were available for alpha diversity.", call. = FALSE)
}

out_root <- file.path(ANALYSES_DIR, "02_Alpha_Diversity")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_csv(alpha, file.path(out_root, "21_frame_alpha_diversity.csv"))

eligibility <- alpha %>%
  count(alpha_model_exclusion_reason, name = "frames") %>%
  mutate(percent = 100 * frames / sum(frames))
write_csv(eligibility, file.path(out_root, "21_alpha_model_eligibility.csv"))

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    observed_common_id_richness, alpha_model_eligible
  ),
  file.path(out_root, "Richness", "21_observed_richness.csv")
)

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    rarefaction_eligible, rarefied_common_id_richness,
    alpha_model_eligible
  ),
  file.path(out_root, "Rarefaction", "21_rarefied_richness.csv")
)

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    shannon_common_id_diversity, alpha_model_eligible
  ),
  file.path(out_root, "Shannon_Wiener", "21_shannon_wiener.csv")
)

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    inverse_simpson_common_id_diversity, alpha_model_eligible
  ),
  file.path(out_root, "Simpson", "21_inverse_simpson.csv")
)

write_csv(
  alpha %>% select(
    filename, dive_id, assigned_biotic_points,
    pielou_evenness, alpha_model_eligible
  ),
  file.path(out_root, "Evenness", "21_pielou_evenness.csv")
)

diversity_summary <- alpha %>%
  summarise(
    point_sampled_frames = n(),
    frames_with_resolved_biota = sum(has_resolved_biota),
    frames_without_resolved_biota = sum(!has_resolved_biota),
    frames_rarefaction_eligible = sum(rarefaction_eligible),
    frames_alpha_model_eligible = sum(alpha_model_eligible),
    median_assigned_biotic_points = median(assigned_biotic_points, na.rm = TRUE),
    median_observed_richness = median(observed_common_id_richness, na.rm = TRUE),
    median_rarefied_richness = median(rarefied_common_id_richness, na.rm = TRUE),
    median_shannon = median(shannon_common_id_diversity, na.rm = TRUE),
    median_inverse_simpson = median(inverse_simpson_common_id_diversity, na.rm = TRUE),
    median_evenness = median(pielou_evenness, na.rm = TRUE)
  )

write_csv(diversity_summary, file.path(out_root, "21_alpha_diversity_summary.csv"))

# Figures intended for inferential interpretation use only frames satisfying
# the predeclared model eligibility criteria. The full descriptive table remains
# available in 21_frame_alpha_diversity.csv.
model_alpha <- alpha %>% filter(alpha_model_eligible)

depth_data <- model_alpha %>% filter(!is.na(depth_m))

if (nrow(depth_data) > 0) {
  p3 <- ggplot(depth_data, aes(x = depth_m, y = observed_common_id_richness)) +
    geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = "Observed biological-group richness across depth",
      subtitle = paste0("Model-eligible frames only; n = ", nrow(depth_data)),
      x = "Depth (m; positive downward)",
      y = "Observed common-ID richness"
    ) +
    theme_bw(base_size = 11)
  save_figure(p3, "Fig_03_Richness_Depth")

  shannon_depth <- depth_data %>% filter(!is.na(shannon_common_id_diversity))
  if (nrow(shannon_depth) > 0) {
    p4 <- ggplot(shannon_depth, aes(x = depth_m, y = shannon_common_id_diversity)) +
      geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = "Shannon diversity across depth",
        subtitle = paste0("Model-eligible frames only; n = ", nrow(shannon_depth)),
        x = "Depth (m; positive downward)",
        y = "Shannon-Wiener diversity (H')"
      ) +
      theme_bw(base_size = 11)
    save_figure(p4, "Fig_04_Shannon_Depth")
  }
}

substrate_data <- model_alpha %>%
  filter(
    !is.na(shannon_common_id_diversity),
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
      subtitle = paste0(
        "Model-eligible frames; substrate classes with at least ",
        MIN_GROUP_N, " frames"
      ),
      x = "Frame substrate class",
      y = "Shannon-Wiener diversity (H')"
    ) +
    theme_bw(base_size = 10)

  save_figure(p5, "Fig_05_Shannon_Substrate", width = 10, height = 9)
}

cat("Point-sampled frames: ", format(nrow(alpha), big.mark = ","), "\n", sep = "")
cat("Frames with resolved biota: ", format(sum(alpha$has_resolved_biota), big.mark = ","), "\n", sep = "")
cat("Rarefaction-eligible frames (>= ", RAREFACTION_N, " Biotic points): ",
    format(sum(alpha$rarefaction_eligible), big.mark = ","), "\n", sep = "")
cat("Alpha-model-eligible frames: ", format(sum(alpha$alpha_model_eligible), big.mark = ","), "\n", sep = "")
cat("Output: ", file.path(out_root, "21_frame_alpha_diversity.csv"), "\n", sep = "")
