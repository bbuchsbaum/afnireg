context("3dDeconvolve nodata tests")

library(fmrireg)
library(fmridesign)
library(fmrihrf)
library(afnireg)

test_that("simple one-regressor model with -nodata works", {
  skip_if_not(file.exists("~/abin/3dDeconvolve"), 
              "3dDeconvolve not found in ~/abin")
  
  # Create minimal event table with one stimulus
  simple_design <- data.frame(
    stim = factor("A"),
    onset = 10,  # Single event at 10 seconds
    run = 1
  )
  
  # Set up sampling frame (100 time points, TR=2)
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Create simple event model with AFNI HRF (needed for -nodata)
  emodel <- event_model(onset ~ afni_hrf(stim, basis = "SPMG1"), 
                       data = simple_design, 
                       block = ~ run, 
                       sampling_frame = sframe)
  
  # Minimal baseline (just intercept/constant)
  bmodel <- baseline_model(basis = "constant", sframe = sframe)
  
  # Create dataset with placeholder files (won't be used with -nodata)
  dset <- fmri_dataset(scans = "dummy.nii", 
                      mask = "dummy_mask.nii",
                      TR = 2,
                      run_length = 100,
                      event_table = simple_design,
                      base_path = ".",
                      dummy_mode = TRUE)
  
  # Combine into fmri_model
  fmodel <- fmri_model(emodel, bmodel, dset)
  
  # Generate AFNI command with -nodata
  alm <- afni_lm(fmodel, dset, nodata = c(100, 2))
  
  # Check that the command contains -nodata
  expect_true(grepl("-nodata 100 2", alm$cmd$cmd))
  
  # Check that no input files are specified
  expect_null(alm$cmd$cmdlines$input)
  expect_null(alm$cmd$cmdlines$mask)
  
  # Test command generation and execution
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "afni_nodata_test")
  
  # Run the command (should create design matrix)
  result <- tryCatch({
    run(alm, outdir = test_dir, execute = TRUE)
    TRUE
  }, error = function(e) {
    message("Error running 3dDeconvolve: ", e$message)
    FALSE
  })
  
  # Check if command executed successfully
  expect_true(result, "3dDeconvolve command should execute without errors")
  
  # Check for expected output files
  if (result && dir.exists(test_dir)) {
    # Check for design matrix file - default bucket name is "stats"
    xmat_file <- file.path(test_dir, "stats.xmat.1D")
    expect_true(file.exists(xmat_file) || 
                file.exists(file.path(test_dir, "X.xmat.1D")) ||
                file.exists(file.path(test_dir, "Decon.xmat.1D")),
                "Design matrix file should be created")
    
    # Clean up
    unlink(test_dir, recursive = TRUE)
  }
})

test_that("nodata with multiple regressors works", {
  skip_if_not(file.exists("~/abin/3dDeconvolve"), 
              "3dDeconvolve not found in ~/abin")
  
  # Create event table with two conditions
  multi_design <- data.frame(
    stim = factor(c("A", "B", "A", "B")),
    onset = c(10, 30, 50, 70),
    run = 1
  )
  
  # Set up sampling frame (100 time points, TR=2) - matches dummy_mode default
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Create event model with AFNI HRF (needed for -nodata)
  emodel <- event_model(onset ~ afni_hrf(stim, basis = "SPMG1"), 
                       data = multi_design, 
                       block = ~ run, 
                       sampling_frame = sframe)
  
  # Baseline with polynomial drift
  bmodel <- baseline_model(basis = "poly", degree = 2, sframe = sframe)
  
  # Create dataset
  dset <- fmri_dataset(scans = "dummy.nii", 
                      mask = "dummy_mask.nii",
                      TR = 2,
                      run_length = 100,
                      event_table = multi_design,
                      base_path = ".",
                      dummy_mode = TRUE)
  
  # Combine into fmri_model
  fmodel <- fmri_model(emodel, bmodel, dset)
  
  # Generate AFNI command with -nodata
  alm <- afni_lm(fmodel, dset, nodata = c(100, 2), 
                options = list(x1D_stop = TRUE))
  
  # Check command structure
  expect_true(grepl("-nodata 100 2", alm$cmd$cmd))
  expect_true(grepl("-x1D_stop", alm$cmd$cmd))
  
  # Test that we have multiple stimulus labels
  expect_equal(length(alm$cmd$afni_stims), 2)  # Two conditions: A and B
})

test_that("nodata with AFNI native HRF works", {
  skip_if_not(file.exists("~/abin/3dDeconvolve"), 
              "3dDeconvolve not found in ~/abin")
  
  # Create simple event table
  design <- data.frame(
    cond = factor("task"),
    onset = c(10, 40, 70),
    run = 1
  )
  
  # Set up sampling frame
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Create event model with AFNI native HRF
  emodel <- event_model(onset ~ afni_hrf(cond, basis = "GAM"), 
                       data = design, 
                       block = ~ run, 
                       sampling_frame = sframe)
  
  # Minimal baseline
  bmodel <- baseline_model(basis = "constant", sframe = sframe)
  
  # Create dataset
  dset <- fmri_dataset(scans = "dummy.nii", 
                      mask = "dummy_mask.nii",
                      TR = 2,
                      run_length = 100,
                      event_table = design,
                      base_path = ".",
                      dummy_mode = TRUE)
  
  # Combine
  fmodel <- fmri_model(emodel, bmodel, dset)
  
  # Generate command
  alm <- afni_lm(fmodel, dset, nodata = c(100, 2))
  
  # Check that -stim_times is used instead of -stim_file for AFNI HRF
  expect_true(grepl("-stim_times", alm$cmd$cmd))
  expect_true(grepl("GAM", alm$cmd$cmd))
})