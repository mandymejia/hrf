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

# Create facet labels
freq_df$c_label <- ifelse(freq_df$c == 0,
                          "No Undershoot (c=0)",
                          "With Undershoot (c=1/6)")

# Gridline positions
a1_vals <- sort(unique(freq_df$a1))
b1_vals <- sort(unique(freq_df$b1))
a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

# ---- Define zero-value tiles for both facets ----
zeros_df <- data.frame(
  a1 = c(
    3, 3, 4, 5,          # original for c = 0
    3, 3, 4,             # original for c = 1/6
    8, 9, 10, 11, 12,    # c = 0, b1 = 0.5
    4,                   # c = 0, b1 = 2
    10, 11, 12,          # c = 0, b1 = 0.75
    12,                  # c = 0, b1 = 1
    10, 11, 12,          # c = 0, b1 = 2
    8, 9, 10, 11, 12,    # c = 1/6, b1 = 0.5
    10, 11, 12,          # c = 1/6, b1 = 0.75
    11,                  # c = 1/6, b1 = 1
    4                    # ✅ new c = 1/6, b1 = 2
  ),
  b1 = c(
    0.5, 0.75, 0.5, 0.5, # original c = 0
    0.5, 0.75, 0.5,      # original c = 1/6
    rep(0.5, 5),         # c = 0, b1 = 0.5
    2,                   # c = 0, b1 = 2
    rep(0.75, 3),        # c = 0, b1 = 0.75
    1,                   # c = 0, b1 = 1
    rep(2, 3),           # c = 0, b1 = 2
    rep(0.5, 5),         # c = 1/6, b1 = 0.5
    rep(0.75, 3),        # c = 1/6, b1 = 0.75
    1,                   # c = 1/6, b1 = 1
    2                    # ✅ new c = 1/6, b1 = 2
  ),
  c_label = c(
    rep("No Undershoot (c=0)", 4),      # original c = 0
    rep("With Undershoot (c=1/6)", 3),  # original c = 1/6
    rep("No Undershoot (c=0)", 5),      # c = 0, b1 = 0.5
    "No Undershoot (c=0)",              # c = 0, b1 = 2
    rep("No Undershoot (c=0)", 3),      # c = 0, b1 = 0.75
    "No Undershoot (c=0)",              # c = 0, b1 = 1
    rep("No Undershoot (c=0)", 3),      # c = 0, b1 = 2
    rep("With Undershoot (c=1/6)", 5),  # c = 1/6, b1 = 0.5
    rep("With Undershoot (c=1/6)", 3),  # c = 1/6, b1 = 0.75
    "With Undershoot (c=1/6)",          # c = 1/6, b1 = 1
    "With Undershoot (c=1/6)"           # ✅ new c = 1/6, b1 = 2
  ),
  rel_freq = 0,
  freq = 0
)
# Merge into the main data (keeps all existing tiles)
freq_df_aug <- bind_rows(freq_df, zeros_df)

# ---- Plot ----
p <- ggplot(freq_df_aug, aes(x = a1, y = b1, fill = rel_freq)) +
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
    'Relative Frequency', option = 'E',
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, max(freq_df$rel_freq))
  ) +
  scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
  scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
  facet_grid(. ~ c_label) +
  theme_few()

print(p)
