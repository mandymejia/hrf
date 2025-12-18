
library(dplyr)
library(ggplot2)
library(ggExtra)

regularize_allHRFs_result <- readRDS("/Volumes/LaCie/root/hrf_adaptation/validation/regularize_allHRFs/regularize_allHRFs_result_tfMRI_MOTOR_LR_1080s.rds")
regularize_allHRFs_result <- readRDS("~/Downloads/regularize_bonferroni_motor_lr_1080s.rds")
regularize_allHRFs_result <- readRDS("~/Downloads/regularize_bh_motor_lr_1080s.rds")

raw_data <- regularize_allHRFs_result[["best_params_df"]] %>%
  select(voxel, subject, a1_raw = a1, mask)

group_means <- regularize_allHRFs_result[["regularized_params"]][["a1"]][["results_OLS"]] %>%
  select(voxel, subject, a1_group_mean = param_mean0, use)

plot_data <- inner_join(raw_data, group_means, by = c("voxel", "subject")) %>%
  filter(use)

ggplot(plot_data, aes(x = a1_group_mean, y = a1_raw)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Group Mean (within subject mask)",
       y = "Subject Raw Estimates",
       title = "Parameter a1: Raw vs Group Mean")


# ============================================

subject_id <- "527"
param <- "a1"


# gt raw estimates and group means for the selected parameter
raw_data <- regularize_allHRFs_result[["best_params_df"]] %>%
  select(voxel, subject, mask, param_raw = !!sym(param))

group_means <- regularize_allHRFs_result[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

# merge and filter to ONE subject within their mask
plot_data <- inner_join(raw_data, group_means, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)

# ask ai to make me the plot for whatever looks good
p <- ggplot(plot_data, aes(x = param_group_mean, y = param_raw)) +
  geom_jitter(width = 0, alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth() +
  geom_smooth(method = "lm", color = "purple") +
  labs(
    x = "Group Mean (within subject mask)",
    y = "Subject Raw Estimates",
    title = paste0("(bh 0.01) Parameter ", param, ": Subject ", subject_id, " (n=", nrow(plot_data), " voxels)")
  ) +
  theme_minimal()

g <- ggMarginal(p, type = "histogram", margins = "y")
ggsave("~/Documents/Github/hrf-z/dev/test_plots/regularize_allHRFs/mean_vs_raw.png", g, width = 7, height = 5, bg="white")
