library(dplyr)
library(ggplot2)
library(ggthemes)

# Use your actual dataset
rounded_params <- regularize_allHRFs_result[["rounded_params"]][["rounded_intercept_slope"]]
# rounded_params <- regularize_allHRFs_result[["best_params_df"]]

hrf_grid <- regularize_allHRFs_result[["hrf_grid"]]


rounded_params <- rounded_params %>%
  mutate(c = round(c, 7))

# Get the complete grid of a1, b1, c combinations from hrf_grid
complete_grid <- hrf_grid %>%
  select(a1, b1, c) %>%
  distinct() %>%
  mutate(c = round(c, 7))  # Round here too



# Aggregate frequencies from rounded_params
freq_df <- rounded_params %>%
  group_by(a1, b1, c) %>%
  summarise(freq = n(), .groups = "drop")

# Left join to get zeros for missing combinations
freq_df <- complete_grid %>%
  left_join(freq_df, by = c("a1", "b1", "c")) %>%
  mutate(freq = ifelse(is.na(freq), 0, freq))

# Compute relative frequency for fill (legend only)
freq_df <- freq_df %>%
  mutate(rel_freq = freq / sum(freq))

# Create facet labels
freq_df$c_label <- ifelse(freq_df$c == 0,
                          "No Undershoot (c=0)",
                          "With Undershoot (c=1/6)")

# Gridline positions - now based on the complete grid
a1_vals <- sort(unique(complete_grid$a1))
b1_vals <- sort(unique(complete_grid$b1))
a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

# ---- Plot ----
p <- ggplot(freq_df, aes(x = a1, y = b1, fill = rel_freq)) +
  geom_hline(yintercept = b1_grid, alpha = 0.1) +
  geom_vline(xintercept = a1_grid, alpha = 0.1) +
  geom_tile(alpha = 0.8) +
  geom_text(aes(label = ifelse(freq == 0, "0", freq))) +
  geom_rect(
    data = subset(freq_df, a1 == 6 & b1 == 1 & abs(c - 1/6) < 1e-10),
    aes(xmin = a1 - 0.5, xmax = a1 + 0.5,
        ymin = b1 - 0.125, ymax = b1 + 0.125),
    color = "black", size = 1.2, fill = NA
  ) +
  scale_fill_viridis_c(
    'Relative \nFrequency', option = 'A',
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max(freq_df$rel_freq))
  ) +
  scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
  scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
  facet_grid(. ~ c_label) +
  theme_few()

print(p)
