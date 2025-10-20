library(dplyr)
library(ggplot2)
library(ggthemes)

# Use your actual dataset
rounded_params <- regularize_allHRFs_result[["rounded_params"]]

# Aggregate frequencies
freq_df <- rounded_params %>%
  group_by(a1, b1, c) %>%
  summarise(freq = n(), .groups = "drop")

# Compute relative frequency for fill (legend only)
freq_df <- freq_df %>%
  mutate(rel_freq = freq / sum(freq))

# Create labels for facets
freq_df$c_label <- ifelse(freq_df$c == 0,
                          "No Undershoot (c=0)",
                          "With Undershoot (c=1/6)")

# Get gridline positions (same logic as in plot_param_grid_metrics)
a1_vals <- sort(unique(freq_df$a1))
b1_vals <- sort(unique(freq_df$b1))
a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

# p <- ggplot(freq_df, aes(x = a1, y = b1, fill = rel_freq)) +
#   geom_hline(yintercept = b1_grid, alpha = 0.1) +
#   geom_vline(xintercept = a1_grid, alpha = 0.1) +
#   geom_tile(alpha = 0.8) +
#   geom_text(aes(label = freq)) +
#   scale_fill_viridis_c('Relative Frequency', option = 'D',
#                        labels = scales::percent_format(accuracy = 1)) +
#   scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
#   scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
#   facet_grid(. ~ c_label) +
#   theme_few()

# Box highlight
p <- ggplot(freq_df, aes(x = a1, y = b1, fill = rel_freq)) +
  geom_hline(yintercept = b1_grid, alpha = 0.1) +
  geom_vline(xintercept = a1_grid, alpha = 0.1) +
  geom_tile(alpha = 0.8) +
  geom_text(aes(label = freq)) +
  geom_rect(
    data = subset(freq_df, a1 == 6 & b1 == 1 & abs(c - 1/6) < 1e-10),
    aes(xmin = a1 - 0.5, xmax = a1 + 0.5,
        ymin = b1 - 0.125, ymax = b1 + 0.125),
    color = "black", size = 1.2, fill = NA
  ) +
  scale_fill_viridis_c('Relative Frequency', option = 'D',
                       labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
  scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
  facet_grid(. ~ c_label) +
  theme_few()



print(p)
