
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

# ============================================ LR AND RL

regularize_result_LR <- readRDS("~/Downloads/regularize_bh_motor_lr_1080s.rds")
regularize_result_RL <- readRDS("~/Downloads/regularize_bh_motor_rl_1085s.rds")

subject_id <- "527"
param <- "a1"

# Extract LR data
raw_data_LR <- regularize_result_LR[["best_params_df"]] %>%
  select(voxel, subject, mask, param_raw = !!sym(param)) %>%
  mutate(run = "LR")

group_means_LR <- regularize_result_LR[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

plot_data_LR <- inner_join(raw_data_LR, group_means_LR, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)

# Extract RL data
raw_data_RL <- regularize_result_RL[["best_params_df"]] %>%
  select(voxel, subject, mask, param_raw = !!sym(param)) %>%
  mutate(run = "RL")

group_means_RL <- regularize_result_RL[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

plot_data_RL <- inner_join(raw_data_RL, group_means_RL, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)

# Combine LR and RL
plot_data <- bind_rows(plot_data_LR, plot_data_RL)

# Create plot with separate linear fits by run
p <- ggplot(plot_data, aes(x = param_group_mean, y = param_raw)) +
  geom_jitter(width = 0, alpha = 0.2, aes(color = run)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", aes(color = run, group = run)) +
  labs(
    x = "Group Mean (within subject mask)",
    y = "Subject Raw Estimates",
    title = paste0("(bh 0.01) Parameter ", param, ": Subject ", subject_id,
                   " (LR: n=", nrow(plot_data_LR), ", RL: n=", nrow(plot_data_RL), " voxels)")
  ) +
  theme_minimal()
ggMarginal(p, type = "histogram", margins = "y", groupFill = TRUE, groupColour = TRUE)
g <- ggMarginal(p, type = "histogram", margins = "y", groupFill = TRUE, groupColour = TRUE)
ggsave("~/Documents/Github/hrf-z/dev/test_plots/regularize_allHRFs/mean_vs_raw_LR_RL.png",
       g, width = 7, height = 5, bg = "white")
# ============================================ LR AND RL (DROP EXTREMES)


regularize_result_LR <- readRDS("~/Downloads/regularize_bonferroni_motor_lr_1080s.rds")
subject_id <- "527"
param <- "a1"

# Get grid extremes - remove 2 steps from each side
hrf_grid <- regularize_result_LR[["hrf_grid"]]  # Same for both
a1_vals <- sort(unique(hrf_grid$a1))
b1_vals <- sort(unique(hrf_grid$b1))

# Exclude first 2 and last 2 values
a1_exclude <- c(a1_vals[1:2], a1_vals[(length(a1_vals)-1):length(a1_vals)])
b1_exclude <- c(b1_vals[1:2], b1_vals[(length(b1_vals)-1):length(b1_vals)])

cat("Will exclude a1 values:", paste(a1_exclude, collapse = ", "), "\n")
cat("Will exclude b1 values:", paste(b1_exclude, collapse = ", "), "\n")

# Extract LR data

raw_data_LR <- regularize_result_LR[["best_params_df"]] %>%
  mutate(param_raw = !!sym(param), run = "LR") %>%
  select(voxel, subject, mask, a1, b1, param_raw, run)

group_means_LR <- regularize_result_LR[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

plot_data_LR <- inner_join(raw_data_LR, group_means_LR, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)
rm(regularize_result_LR); gc()
# Extract RL data
regularize_result_RL <- readRDS("~/Downloads/regularize_bonferroni_motor_rl_1085s.rds")
raw_data_RL <- regularize_result_RL[["best_params_df"]] %>%
  mutate(param_raw = !!sym(param), run = "RL") %>%
  select(voxel, subject, mask, a1, b1, param_raw, run)

group_means_RL <- regularize_result_RL[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

plot_data_RL <- inner_join(raw_data_RL, group_means_RL, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)

# Combine LR and RL
plot_data <- bind_rows(plot_data_LR, plot_data_RL)

cat("Before filtering: LR =", nrow(plot_data_LR), ", RL =", nrow(plot_data_RL),
    ", Total =", nrow(plot_data), "voxels\n")

# WIPING OUT: Remove voxels at grid extremes (2 steps from each side)
plot_data_filtered <- plot_data %>%
  filter(!(a1 %in% a1_exclude)) %>%
  filter(!(b1 %in% b1_exclude))

cat("After filtering:", nrow(plot_data_filtered), "voxels\n")
cat("Removed:", nrow(plot_data) - nrow(plot_data_filtered), "voxels\n")

# Create plot with LR and RL in different colors
p <- ggplot(plot_data_filtered, aes(x = param_group_mean, y = param_raw)) +
  geom_jitter(width = 0, alpha = 0.2, aes(color = run)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", aes(color = run, group = run)) +
  labs(
    x = "Group Mean (within subject mask)",
    y = "Subject Raw Estimates",
    title = paste0("Parameter ", param, ": Subject ", subject_id,
                   " (2 Steps Removed, n=", nrow(plot_data_filtered), " voxels)")
  ) +
  theme_minimal()
ggMarginal(p, type = "histogram", margins = "y", groupFill = TRUE, groupColour = TRUE)
g <- ggMarginal(p, type = "histogram", margins = "y", groupFill = TRUE, groupColour = TRUE)
ggsave("~/Documents/Github/hrf-z/dev/test_plots/regularize_allHRFs/mean_vs_raw_LR_RL_filtered_2steps.png",
       g, width = 7, height = 5, bg = "white")
# ============================================ LR  or RL (DROP EXTREMES)

# Load your result
regularize_result <- readRDS("~/Downloads/regularize_bh_motor_lr_1080s.rds")

subject_id <- "527"
param <- "a1"

# Get grid extremes - remove 2 steps from each side
hrf_grid <- regularize_result[["hrf_grid"]]
a1_vals <- sort(unique(hrf_grid$a1))
b1_vals <- sort(unique(hrf_grid$b1))

# Exclude first 2 and last 2 values
a1_exclude <- c(a1_vals[1:2], a1_vals[(length(a1_vals)-1):length(a1_vals)])
b1_exclude <- c(b1_vals[1:2], b1_vals[(length(b1_vals)-1):length(b1_vals)])

cat("Will exclude a1 values:", paste(a1_exclude, collapse = ", "), "\n")
cat("Will exclude b1 values:", paste(b1_exclude, collapse = ", "), "\n")

# Get raw estimates and group means
raw_data <- regularize_result[["best_params_df"]] %>%
  mutate(param_raw = !!sym(param)) %>%
  select(voxel, subject, mask, a1, b1, param_raw)

group_means <- regularize_result[["regularized_params"]][[param]][["results_OLS"]] %>%
  select(voxel, subject, use, param_group_mean = param_mean0)

# Merge and filter to ONE subject within their mask
plot_data <- inner_join(raw_data, group_means, by = c("voxel", "subject")) %>%
  filter(subject == subject_id, use)

cat("Before filtering:", nrow(plot_data), "voxels\n")

# WIPING OUT: Remove voxels at grid extremes (2 steps from each side)
plot_data_filtered <- plot_data %>%
  filter(!(a1 %in% a1_exclude)) %>%
  filter(!(b1 %in% b1_exclude))

cat("After filtering:", nrow(plot_data_filtered), "voxels\n")
cat("Removed:", nrow(plot_data) - nrow(plot_data_filtered), "voxels\n")

# Create plot
p <- ggplot(plot_data_filtered, aes(x = param_group_mean, y = param_raw)) +
  geom_jitter(width = 0, alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", color = "purple") +
  labs(
    x = "Group Mean (within subject mask)",
    y = "Subject Raw Estimates",
    title = paste0("Parameter ", param, ": Subject ", subject_id,
                   " (2 Steps Removed, n=", nrow(plot_data_filtered), " voxels)")
  ) +
  theme_minimal()

g <- ggMarginal(p, type = "histogram", margins = "y")
ggsave("~/Documents/Github/hrf-z/dev/test_plots/regularize_allHRFs/mean_vs_raw_filtered_2steps.png",
       g, width = 7, height = 5, bg = "white")

# ====================



subject_raw <- regularize_result[["best_params_df"]] %>%
  filter(subject == "12") %>%
  select(voxel, a1)

xii_template <- attr(regularize_result, "xii")


full_vector <- rep(NA, nrow(as.matrix(xii_template)))
full_vector[subject_raw$voxel] <- subject_raw$a1


raw_a1_xii <-ciftiTools::newdata_xifti(xii_template, full_vector)
raw_a1_xii<-smooth_xifti(raw_a1_xii, surf_FWHM=8)

plot(raw_a1_xii, zlim=c(3,12))
