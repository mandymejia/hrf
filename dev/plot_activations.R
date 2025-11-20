library(here)
library(ciftiTools)
devtools::load_all("~/Documents/Github/hrf-z")


ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')


session_data <- readRDS(here("dev", "fixtures", "session_data_4s", paste0("session_data_", "motor", "_lr_4s.rds")))
BOLD_xii <- ciftiTools::read_cifti(session_data[["BOLD_files"]][[1]],
                                   brainstructures = c("left", "right"),
                                   resamp_res = 10000)
# Fit_workingHRF

design_2d <- make_design(
  EVs = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR = 0.72,
  onset = TRUE,
  offset = TRUE
)

plot(design_2d)

design_array <- array(design_2d$design, dim = c(dim(design_2d$design), 2))
design_array[,,2] <- design_array[,,2] + rnorm(length(design_array[,,2]))

glm_result_1 <- multiGLM(
  BOLD = BOLD_xii,
  brainstructures = c("left", "right"),
  design = design_array,
  design_canonical = design_array[,,1],
  nuisance = as.matrix(read.table(session_data[["nuisance_files"]][[1]], header=FALSE)),
  TR = 0.72,
  hpf = 0.01,
  resamp_res = NULL
)

plot(glm_result_1$Fstat_xii)
plot(glm_result_1$pvalF_xii)

# Fit_ALLHRFS


#### LETS SCALE THE BOLD FILE
BOLD_xii <- ciftiTools::read_cifti(session_data[["BOLD_files"]][[1]],
                                   brainstructures = c("left", "right"),
                                   resamp_res = 10000)
# template_xii <- ciftiTools::convert_xifti(BOLD_xii, to = "dscalar") # Used with and without template
mat <- as.matrix(BOLD_xii)

## scale
mu <- rowMeans(mat)
mat_psc <- ((mat - mu) / mu) * 100

scaled_BOLD_xii <- ciftiTools::newdata_xifti(BOLD_xii, mat_psc) # template_xii or BOLD_xii tried here
BOLD_xii <- scaled_BOLD_xii

small_hrf_grid <- default_hrf_grid#[1:3, ]
d1 <- make_design(
  EVs = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR = 0.72,
  dHRF = 0,
  onset = TRUE,
  offset = TRUE,
  taper_start = NULL,
  a1 = small_hrf_grid$a1[1],
  b1 = small_hrf_grid$b1[1],
  c  = small_hrf_grid$c[1],
  a2 = small_hrf_grid$a2[1],
  b2 = small_hrf_grid$b2[1]
)

d2 <- make_design(
  EVs = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR = 0.72,
  dHRF = 0,
  onset = TRUE,
  offset = TRUE,
  taper_start = NULL,
  a1 = small_hrf_grid$a1[2],
  b1 = small_hrf_grid$b1[2],
  c  = small_hrf_grid$c[2],
  a2 = small_hrf_grid$a2[2],
  b2 = small_hrf_grid$b2[2]
)

d3 <- make_design(
  EVs = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR = 0.72,
  dHRF = 0,
  onset = TRUE,
  offset = TRUE,
  taper_start = NULL,
  a1 = small_hrf_grid$a1[3],
  b1 = small_hrf_grid$b1[3],
  c  = small_hrf_grid$c[3],
  a2 = small_hrf_grid$a2[3],
  b2 = small_hrf_grid$b2[3]
)


d4 <- make_design(
  EVs = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR = 0.72,
  dHRF = 0,
  onset = TRUE,
  offset = TRUE,
  taper_start = NULL,
  a1 = small_hrf_grid$a1[24],
  b1 = small_hrf_grid$b1[24],
  c  = small_hrf_grid$c[24],
  a2 = small_hrf_grid$a2[24],
  b2 = small_hrf_grid$b2[24]
)

design_3D <- array(NA, dim = c(ncol(BOLD_xii), ncol(d1$design), 4))
design_3D[,,1] <- d1$design
design_3D[,,2] <- d2$design
design_3D[,,3] <- d3$design
design_3D[,,4] <- d4$design

devtools::load_all("~/Documents/Github/hrf-z")


glm_result <- multiGLM(
  BOLD = BOLD_xii,
  brainstructures = c("left", "right"),
  design = design_3D,
  scrub = NULL,
  nuisance = as.matrix(read.table(session_data[["nuisance_files"]][[1]], header = FALSE)),
  TR = 0.72,
  hpf = 0.01,
  resamp_res = NULL
)

plot(glm_result[["GLMs"]][[4]][["betas"]], idx = 1, zlim = c(-1, 1))
plot(glm_result[["GLMs"]][[4]][["betas"]], idx = 6, zlim = c(-1, 1))
plot(glm_result[["GLMs"]][[4]][["betas"]], idx = 1, zlim = c(-7.37, 10.8))
plot(glm_result[["GLMs"]][[4]][["betas"]], idx = 6, zlim = c(-1350,904))
################################################################################
# CLASSICAL GLM
library(BayesfMRI)
################################################################################

# shared: build design + nuisance
design <- BayesfMRI::make_design(
  EVs   = session_data[["EVs_list"]][[1]],
  nTime = ncol(BOLD_xii),
  TR    = 0.72
)

# 1) Classical GLM (massive univariate), minimal
bglm_classic <- BayesGLM(
  BOLD          = BOLD_xii,
  design        = design$design,     # 2D matrix (time x fields)
  nuisance      = as.matrix(read.table(session_data[["nuisance_files"]][[1]], header = FALSE)),
  TR            = 0.72,
  hpf           = 0.01 ,
  brainstructures = c("left", "right"),
  Bayes         = FALSE,             # <- classical
  verbose       = 0
)

# final “picture” = activations map; use FWER or FDR for classical
act_classic <- activations(bglm_classic, Bayes = FALSE, alpha = 0.05, correction = "FDR", verbose = 0)
plot(act_classic, idx = design$field_names[1], what = "surface", title = "Classical GLM activations")

# stat map instead of thresholded activations:
plot(bglm_classic, Bayes = FALSE, idx = design$field_names[1], what = "surface", title = "Classical GLM beta")

# 2) Spatial Bayesian GLM (looks cleaner), minimal
bglm_bayes <- BayesGLM(
  BOLD            = BOLD_xii,
  design          = design$design,
  nuisance        = as.matrix(read.table(session_data[["nuisance_files"]][[1]], header = FALSE)),
  TR              = 0.72,
  hpf             = 0.01,
  brainstructures = c("left", "right"),
  Bayes           = TRUE,    # <- spatial model
  ar_order        = 3,       # prewhitening; set 0 to skip
  ar_smooth       = 0,       # keep 0 for bare-bones
  verbose         = 0
)

# thresholded activations (Bayesian excursion sets)
act_bayes <- activations(bglm_bayes, Bayes = TRUE, alpha = 0.05, gamma = 0, verbose = 0)
plot(act_bayes, idx = design$field_names[1], what = "surface", title = "Spatial Bayes activations")

# subcortex slice
plot(act_bayes, idx = design$field_names[1], what = "volume", slices = 27:30, n_slices = 4)
