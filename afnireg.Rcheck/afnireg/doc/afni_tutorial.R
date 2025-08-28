## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  message = FALSE,
  warning = FALSE,
  error = TRUE # Allow errors to be displayed for AFNI-specific checks
)
library(tibble)
library(fmrireg)
# Check if AFNI is available - needed for AFNI-specific HRFs
# You might replace this with a more robust check, e.g., system("which 3dDeconvolve")
afni_available <- FALSE # Set to TRUE if AFNI is installed and in PATH

## ----setup_simple_afni--------------------------------------------------------
TR <- 2
cond <- c("face", "scene", "tool", "object")
NSTIM <- length(cond) * 4

set.seed(123)
simple_design <- data.frame(
  stim = factor(sample(rep(cond, 4))),
  ISI = sample(10:20, NSTIM, replace = TRUE), # Using wider ISI from previous fixes
  run = rep(1, NSTIM),
  trial = factor(1:NSTIM)
)

simple_design$onset <- cumsum(c(0, simple_design$ISI[-NSTIM] + 2))
sframe <- sampling_frame(blocklens = 140, TR = TR)

# Define a contrast
con1 <- pair_contrast(~ stim == "face", ~stim == "scene", name = "face_scene")

# Create the event model
emodel <- event_model(onset ~ hrf(stim, contrasts = con1), 
                     data = simple_design, 
                     block = ~ run, 
                     sampling_frame = sframe)

# Create the baseline model (B-spline drift)
bmodel <- baseline_model(basis = "bs", degree = 5, sframe = sframe)

# Combine into a full fmri_model
fmodel <- fmri_model(emodel, bmodel)

# Define the dataset (using placeholder file names)
dset <- fmri_dataset(scans = "scan01.nii.gz",
                     mask = "mask.nii.gz",
                     TR = TR,
                     run_length = 140,
                     event_table = simple_design,
                     base_path = ".") # Assuming files are in the working directory

## ----create_afni_spec---------------------------------------------------------
alm_spec <- afni_lm(fmodel, dset)
print(alm_spec)

## ----afni_polort_example, eval=FALSE------------------------------------------
# # Example: No drift in bmodel, use AFNI polort
# bmodel_no_drift <- baseline_model(basis = "constant", sframe = sframe) # Only intercept
# fmodel_no_drift <- fmri_model(emodel, bmodel_no_drift)
# 
# # Specify polort=2 for quadratic drift
# alm_spec_polort <- afni_lm(fmodel_no_drift, dset, polort = 2)
# print(alm_spec_polort)

## ----afni_options_example-----------------------------------------------------
# Example: Request statistic and coefficient buckets, disable t-stat output
alm_spec_opts <- afni_lm(fmodel, dset, 
                         options = list(bucket = "stats_afni", 
                                        cbucket = "coefs_afni",
                                        tout = FALSE))
print(alm_spec_opts)

## ----afni_native_hrf, eval=afni_available-------------------------------------
# # Example using AFNI's BLOCK function for convolution
# # Note: Using afni_hrf requires AFNI to be installed for 3dDeconvolve to run
# emodel_afni_block <- event_model(
#   onset ~ afni_hrf(stim, basis = "BLOCK(2,1)"), # BLOCK(duration, penalty)
#   data = simple_design,
#   block = ~ run,
#   sampling_frame = sframe
# )
# 
# fmodel_afni <- fmri_model(emodel_afni_block, bmodel)
# alm_spec_afni_native <- afni_lm(fmodel_afni, dset)
# print(alm_spec_afni_native)
# 
# # Check the command: It should now use -stim_times instead of -stim_file
# # and specify 'BLOCK(2,1)' as the basis function.

## ----afni_native_hrf_msg, eval=!afni_available, echo=FALSE--------------------
cat("Skipping AFNI-native HRF example as AFNI availability check is FALSE.")

## ----afni_censor_example------------------------------------------------------
# Example: Censor the first 5 scans
censor_vec <- rep(1, 140) # Start with all scans included (1)
censor_vec[1:5] <- 0       # Censor first 5 scans (0)

alm_spec_censored <- afni_lm(fmodel, dset, censor = censor_vec)
print(alm_spec_censored)

## ----afni_run, eval=FALSE-----------------------------------------------------
# # This command will:
# # 1. Create the output directory 'glm_afni_output'.
# # 2. Change into that directory.
# # 3. Write all necessary files:
# #    - Stimulus files (.1D) for each regressor if using standard hrf()
# #    - Stimulus timing files (.1D) if using afni_hrf()
# #    - Nuisance regressor files (.1D)
# #    - GLT files (.txt)
# #    - Censor file (censor.1D) if specified
# #    - The 3ddeconvolve.sh script containing the full command
# # 4. Execute the 3ddeconvolve.sh script.
# # 5. Change back to the original directory.
# 
# run(alm_spec_opts, outdir = "glm_afni_output")
# 
# # To only generate the script and files without running:
# # run(alm_spec_opts, outdir = "glm_afni_output", execute = FALSE)

