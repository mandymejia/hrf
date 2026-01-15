library(ggplot2)

result_LR <- readRDS("/Volumes/LaCie/root/hrf_adaptation/validation_bonferroni_0.05/fit_workingHRF/fit_workingHRF_result_tfMRI_MOTOR_LR_1080s.rds")
result_RL <- readRDS("/Volumes/LaCie/root/hrf_adaptation/validation_bonferroni_0.05/fit_workingHRF/fit_workingHRF_result_tfMRI_MOTOR_RL_1085s.rds")

n_active_LR <- colSums(result_LR$activation_masks$masks, na.rm = TRUE)
n_active_RL <- colSums(result_RL$activation_masks$masks, na.rm = TRUE)

# Create data frame
plot_data <- data.frame(
  n_active = c(n_active_LR, n_active_RL),
  run = rep(c("LR", "RL"), c(length(n_active_LR), length(n_active_RL)))
)

# Overlapping histogram
ggplot(plot_data, aes(x = n_active, fill = run)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30, color = "black") +
  scale_fill_manual(values = c("LR" = "#2c7fb8", "RL" = "#d95f02")) +
  labs(
    x = "Number of Active Locations per Subject",
    y = "Number of Subjects",
    title = "Motor Task: Active Locations per Subject (LR vs RL)",
    fill = "Run"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Print summary stats
cat("LR: mean =", mean(n_active_LR), ", sd =", sd(n_active_LR), "\n")
cat("RL: mean =", mean(n_active_RL), ", sd =", sd(n_active_RL), "\n")
