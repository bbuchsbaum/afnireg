#!/usr/bin/env Rscript

# Test mixed model with both R-based hrf() and AFNI-native afni_hrf() terms

library(fmridesign)
library(fmrihrf)
library(afnireg)

# Create a simple design with two conditions and a covariate
set.seed(123)
design <- data.frame(
  condition = factor(rep(c("A", "B"), each = 10)),
  RT = rnorm(20, mean = 0.5, sd = 0.2),
  run = rep(1:2, each = 10)
)

# Generate onsets
design$onset <- seq(1, 200, length.out = 20)

# Set up sampling frame
sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)

# Create baseline model
bmodel <- baseline_model(basis = "bs", degree = 3, sframe = sframe)

# Define contrasts for AFNI term
con_A_vs_B <- contrast(~ condition.A - condition.B, name = "A_vs_B")

# Test 1: Mixed model with both hrf() and afni_hrf()
cat("\n=== TEST 1: Mixed Model (hrf + afni_hrf) ===\n")
emodel_mixed <- event_model(
  onset ~ hrf(RT, basis = "spmg1") + afni_hrf(condition, basis = "SPMG1", contrasts = list(con_A_vs_B)),
  data = design,
  block = ~ run,
  sampling_frame = sframe
)

# Create dummy dataset
dataset <- fmri_dataset(
  scans = "dummy.nii",
  mask = "dummy_mask.nii", 
  TR = 2,
  run_length = 100,
  event_table = design,
  base_path = ".",
  dummy_mode = TRUE
)

# Create full model
model_mixed <- fmri_model(emodel_mixed, bmodel, dataset)

# Build AFNI command
cat("\nBuilding AFNI command for mixed model...\n")
cmd_mixed <- tryCatch({
  afnireg:::build_decon_command(
    model = model_mixed,
    dataset = list(scans = NULL, mask_file = NULL),
    working_dir = getwd(),
    opts = list(
      nodata = c(100, 2),
      polort = -1,
      iresp = FALSE,
      TR_times = 2,
      censor = NULL,
      x1D_stop = TRUE,
      nofullf_atall = TRUE
    )
  )
}, error = function(e) {
  cat("Error building command:", e$message, "\n")
  NULL
})

if (!is.null(cmd_mixed)) {
  cat("\n--- Mixed Model Results ---\n")
  cat("Number of AFNI stims:", length(cmd_mixed$afni_stims), "\n")
  cat("Number of GLTs:", length(cmd_mixed$glts), "\n")
  cat("GLT names:", paste(cmd_mixed$gltnames, collapse = ", "), "\n")
  
  # Check for R-based terms (should use -stim_file)
  has_stim_file <- grepl("-stim_file", cmd_mixed$cmd)
  cat("Has -stim_file (R-based):", has_stim_file, "\n")
  
  # Check for AFNI terms (should use -stim_times)
  has_stim_times <- grepl("-stim_times", cmd_mixed$cmd)
  cat("Has -stim_times (AFNI-native):", has_stim_times, "\n")
}

# Test 2: Pure R-based model with hrf()
cat("\n\n=== TEST 2: Pure R Model (hrf only) ===\n")
emodel_r <- event_model(
  onset ~ hrf(condition) + hrf(RT, basis = "spmg1"),
  data = design,
  block = ~ run,
  sampling_frame = sframe
)

model_r <- fmri_model(emodel_r, bmodel, dataset)

cmd_r <- tryCatch({
  afnireg:::build_decon_command(
    model = model_r,
    dataset = list(scans = NULL, mask_file = NULL),
    working_dir = getwd(),
    opts = list(
      nodata = c(100, 2),
      polort = -1,
      iresp = FALSE,
      TR_times = 2,
      censor = NULL,
      x1D_stop = TRUE,
      nofullf_atall = TRUE
    )
  )
}, error = function(e) {
  cat("Error building command:", e$message, "\n")
  NULL
})

if (!is.null(cmd_r)) {
  cat("\n--- R Model Results ---\n")
  cat("Number of AFNI stims:", length(cmd_r$afni_stims), "\n")
  cat("Has -stim_file:", grepl("-stim_file", cmd_r$cmd), "\n")
  cat("Has -stim_times:", grepl("-stim_times", cmd_r$cmd), "\n")
}

# Test 3: Pure AFNI model with afni_hrf()
cat("\n\n=== TEST 3: Pure AFNI Model (afni_hrf only) ===\n")
emodel_afni <- event_model(
  onset ~ afni_hrf(condition, basis = "SPMG1", contrasts = list(con_A_vs_B)) + 
          afni_hrf(RT, basis = "SPMG1"),
  data = design,
  block = ~ run,
  sampling_frame = sframe
)

model_afni <- fmri_model(emodel_afni, bmodel, dataset)

cmd_afni <- tryCatch({
  afnireg:::build_decon_command(
    model = model_afni,
    dataset = list(scans = NULL, mask_file = NULL),
    working_dir = getwd(),
    opts = list(
      nodata = c(100, 2),
      polort = -1,
      iresp = FALSE,
      TR_times = 2,
      censor = NULL,
      x1D_stop = TRUE,
      nofullf_atall = TRUE
    )
  )
}, error = function(e) {
  cat("Error building command:", e$message, "\n")
  NULL
})

if (!is.null(cmd_afni)) {
  cat("\n--- AFNI Model Results ---\n")
  cat("Number of AFNI stims:", length(cmd_afni$afni_stims), "\n")
  cat("Number of GLTs:", length(cmd_afni$glts), "\n")
  cat("GLT names:", paste(cmd_afni$gltnames, collapse = ", "), "\n")
  cat("Has -stim_file:", grepl("-stim_file", cmd_afni$cmd), "\n")
  cat("Has -stim_times:", grepl("-stim_times", cmd_afni$cmd), "\n")
}

cat("\n\n=== TEST SUMMARY ===\n")
cat("Mixed model (hrf + afni_hrf): ", 
    if(!is.null(cmd_mixed)) "SUCCESS" else "FAILED", "\n")
cat("R-only model (hrf): ", 
    if(!is.null(cmd_r)) "SUCCESS" else "FAILED", "\n")
cat("AFNI-only model (afni_hrf): ", 
    if(!is.null(cmd_afni)) "SUCCESS" else "FAILED", "\n")