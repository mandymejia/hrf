library(here)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tictoc)
library(ciftiTools)

devtools::load_all("~/Documents/Github/hrf-z")
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
session_data <- readRDS(here(
  "dev", "fixtures", "session_data_4s",
  "session_data_motor_lr_4s.rds"
))



plot_design <- function(x) {
  design_dHRFs <- x
  volume <- 1:nrow(design_dHRFs)

  # Define which HRF types to include and their suffixes
  hrf_types <- list(
    main = "",           # main regressors: no suffix
    dhrf = "_dHRF",      # dhrf regressors: ends with _dHRF
    ddhrf = "_ddHRF"
  )

  # Collect long-format rows for all HRF types
  df_list <- list()

  for (hrf_name in names(hrf_types)) {
    suffix <- hrf_types[[hrf_name]]

    # Get columns that match the suffix ("" means exclude _dHRF/_ddHRF)
    if (suffix == "") {
      regs <- colnames(design_dHRFs)[!grepl("_dHRF|_ddHRF", colnames(design_dHRFs))]
    } else {
      pattern <- paste0(suffix, "$")
      regs <- grep(pattern, colnames(design_dHRFs), value = TRUE)
    }

    df_list[[hrf_name]] <- do.call(rbind, lapply(regs, function(reg) {
      base_task <- sub(paste0(suffix, "$"), "", reg)  # Strip suffix if present
      data.frame(
        Volume = volume,
        Task = base_task,
        HRF = hrf_name,
        Value = design_dHRFs[, reg]
      )
    }))
  }

  df_long <- do.call(rbind, df_list)
  df_long$Task <- factor(df_long$Task, levels = unique(df_long$Task))

  plot_main <- ggplot(data = df_long[df_long$HRF == "main",]) +
    geom_line(aes(x = Volume, y = Value, color = Task), linetype = "solid") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = "Main HRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    facet_grid(Task ~ .)

  plot_dhrf <- ggplot(data = df_long[df_long$HRF == "dhrf",]) +
    geom_line(aes(x = Volume, y = Value, color = Task), linetype = "dashed") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = "dHRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    facet_grid(Task ~ .)

  plot_ddhrf <- ggplot(data = df_long[df_long$HRF == "ddhrf",]) +
    geom_line(aes(x = Volume, y = Value, color = Task), linetype = "twodash") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = "ddHRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    facet_grid(Task ~ .)

  if ("main" %in% names(df_list)) print(plot_main)
  if ("dhrf" %in% names(df_list)) print(plot_dhrf)
  if ("ddhrf" %in% names(df_list)) print(plot_ddhrf)

  invisible(Filter(Negate(is.null), list(
    main = if ("main" %in% names(df_list)) plot_main else NULL,
    dhrf = if ("dhrf" %in% names(df_list)) plot_dhrf else NULL,
    ddhrf = if ("ddhrf" %in% names(df_list)) plot_ddhrf else NULL
  )))
}

EVs <- session_data[["EVs_list"]][[1]]

design_dHRFs <- make_design(
  EVs=EVs, nTime = 284, TR = 0.72, dHRF = 2, onset = TRUE, offset = TRUE
)$design


plot_design(design_dHRFs)
plots <- plot_design(design_dHRFs)
