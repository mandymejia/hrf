library(dplyr)
library(ggplot2)
devtools::load_all("~/Documents/Github/hrf-z")

# Load default HRF grid
data("default_hrf_grid")

# Compute metrics (time_to_peak, FWHM, etc.)
metrics <- compute_hrf_metrics(default_hrf_grid)
hrf_params <- metrics$hrf_params

# Get unique combos with metrics
unique_combos <- hrf_params %>%
  select(a1, b1, c, time_to_peak, FWHM) %>%
  distinct() %>%
  mutate(
    c_label = ifelse(c == 0, "No Undershoot (c=0)", "With Undershoot (c=1/6)")
  )

cat("Unique combinations:", nrow(unique_combos), "\n")
cat("Time-to-peak range:", range(unique_combos$time_to_peak), "\n")
cat("FWHM range:", range(unique_combos$FWHM), "\n\n")

# Print table
print(unique_combos)

# Plot: time_to_peak vs FWHM, faceted by c
p <- ggplot(unique_combos, aes(x = time_to_peak, y = FWHM)) +
  geom_point(aes(color = factor(a1), size = b1), alpha = 0.8) +
  facet_wrap(~c_label) +
  scale_color_viridis_d(name = "a1") +
  scale_size_continuous(name = "b1", range = c(3, 12), breaks = seq(0.5, 2, 0.25)) +
  labs(
    title = "HRF Grid: Time-to-Peak vs FWHM",
    subtitle = paste0("All ", nrow(unique_combos), " unique (a1, b1, c) combinations"),
    x = "Time to Peak (seconds)",
    y = "FWHM (seconds)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(p)

# Plot: time_to_peak vs FWHM, single plot with shape for c
q <- ggplot(unique_combos, aes(x = time_to_peak, y = FWHM)) +
  geom_point(aes(color = factor(a1), size = b1, shape = c_label), alpha = 0.8) +
  scale_color_viridis_d(name = "a1") +
  scale_size_continuous(name = "b1", range = c(3, 12), breaks = seq(0.5, 2, 0.25)) +
  scale_shape_manual(name = "c", values = c("No Undershoot (c=0)" = 16, "With Undershoot (c=1/6)" = 17)) +
  labs(
    title = "HRF Grid: Time-to-Peak vs FWHM",
    subtitle = paste0("All ", nrow(unique_combos), " unique (a1, b1, c) combinations"),
    x = "Time to Peak (seconds)",
    y = "FWHM (seconds)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

print(q)

# Save plots
plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
ggsave(file.path(plot_dir, "hrf_grid_time_to_peak_vs_fwhm.png"), p, width = 10, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "hrf_grid_time_to_peak_vs_fwhm_combined.png"), q, width = 8, height = 6, dpi = 300)

cat("\nPlots saved to:", plot_dir, "\n")
