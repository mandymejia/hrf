library(here)
library(qs)
library(ggplot2)
library(ggthemes)
library(gridExtra)


ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

hrf_grid <- generate_hrf_grid()
all_hrfs_grid <- generate_hrf_grid(peak2_time_max = 50)


a <- plot(all_hrfs_grid)
ggsave("dev/test_plots/hrf_grid/HRFs_all.pdf", a, width = 8, height = 10.5)
b <- plot(hrf_grid)
ggsave("dev/test_plots/hrf_grid/HRFs.pdf", b, width = 8, height = 10.5)

e <-plot(all_hrfs_grid, type = "hrfs_tapered")
ggsave("dev/test_plots/hrf_grid/HRFs_tapered_all.pdf", e, width = 8, height = 10.5)
f <- plot(hrf_grid, type = "hrfs_tapered")
ggsave("dev/test_plots/hrf_grid/HRFs_tapered.pdf", f, width = 8, height = 10.5)

c <- plot(all_hrfs_grid, type = "param_grid")
ggsave("dev/test_plots/hrf_grid/param_grid_all.pdf", c, width = 10, height = 9)
d <- plot(hrf_grid, type = "param_grid")
ggsave("dev/test_plots/hrf_grid/param_grid.pdf", d, width = 10, height = 9)


g <- plot(hrf_grid, type = "single_hrf", hrf_idx = 24, tapered = FALSE)
ggsave("dev/test_plots/hrf_grid/HRF_single_40.pdf", g, width = 6, height = 4)
h <- plot(hrf_grid, type = "single_hrf", hrf_idx = 24, tapered = TRUE)
ggsave("dev/test_plots/hrf_grid/HRF_single_40_tapered.pdf", h, width = 6, height = 4)

