context("AFNI 3dDeconvolve command generation")

test_that("simple single stimulus generates correct command", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe,
    baseline_degree = 0  # Just constant
  )
  
  cmd <- generate_afni_command(model, nodata = c(100, 2))
  
  # Check basic command structure
  expect_has_nodata(cmd$cmd, 100, 2)
  expect_true(grepl("-polort -1", cmd$cmd))
  expect_equal(get_stim_count(cmd$cmd), 1)
  expect_has_stim_times(cmd$cmd, 1)
  expect_true(grepl("SPMG1", cmd$cmd))
})

test_that("two conditions generate two stim_times", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  expect_equal(get_stim_count(cmd$cmd), 2)
  expect_has_stim_times(cmd$cmd, 2)
  expect_equal(length(cmd$afni_stims), 2)
  
  # Check stim labels
  expect_true(grepl("-stim_label 1", cmd$cmd))
  expect_true(grepl("-stim_label 2", cmd$cmd))
})

test_that("multiple runs generate concat option", {
  design <- make_simple_design(n_cond = 2, n_runs = 3)
  sframe <- make_sampling_frame(n_runs = 3)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model, nodata = c(300, 2))
  
  # With multiple runs and nodata, should have concat
  expect_true(grepl("-concat", cmd$cmd))
  expect_true(grepl("1D: 0 100 200", cmd$cmd))  # Run starts at 0, 100, 200
})

test_that("different HRF models generate correct basis strings", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Test different bases
  bases <- list(
    "SPMG1" = "SPMG1\\(1\\)",
    "GAM" = "GAM\\(8.6,0.547\\)",
    "BLOCK" = "BLOCK\\(1,1\\)",
    "SIN" = "SIN\\(1\\)"
  )
  
  for (basis_name in names(bases)) {
    model <- build_test_afni_model(
      design,
      as.formula(sprintf("onset ~ afni_hrf(condition, basis = '%s')", basis_name)),
      sframe
    )
    
    cmd <- generate_afni_command(model)
    expected_pattern <- bases[[basis_name]]
    
    expect_true(grepl(expected_pattern, cmd$cmd),
                info = paste("Failed for basis:", basis_name))
  }
})

test_that("events with durations work correctly", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  design$duration <- 2.5  # Add duration column
  sframe <- make_sampling_frame(n_runs = 1)
  
  # BLOCK basis typically used with durations
  # Note: afni_hrf requires a constant duration, not a variable
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "BLOCK", durations = 2.5),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  expect_true(grepl("BLOCK", cmd$cmd))
  # Check that stim_times files will include durations
  # (actual file content would be checked in stimulus_files tests)
})

test_that("factorial design generates correct number of conditions", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # 2x2 factorial = 4 conditions
  expect_equal(get_stim_count(cmd$cmd), 4)
  expect_has_stim_times(cmd$cmd, 4)
  expect_equal(length(cmd$afni_stims), 4)
})

test_that("baseline model adds ortvec options", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe,
    baseline_degree = 3  # Polynomial baseline
  )
  
  cmd <- generate_afni_command(model)
  
  # Should have ortvec for baseline regressors
  expect_true(grepl("-ortvec", cmd$cmd))
  expect_gt(length(cmd$afni_baseline_mats), 0)
})

test_that("x1D_stop option works", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model, options = list(x1D_stop = TRUE))
  
  expect_true(grepl("-x1D_stop", cmd$cmd))
})

test_that("TR_times option affects timing", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1, TR = 2.5)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model, nodata = c(100, 2.5))
  
  expect_true(grepl("-TR_times 2.5", cmd$cmd))
})

test_that("command execution with 3dDeconvolve works", {
  skip_on_cran()
  skip_if_not(has_3ddeconvolve(), "3dDeconvolve not found")
  
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  # Create afni_lm object
  dset <- model$dataset
  fmodel <- model
  alm <- afni_lm(fmodel, dset, nodata = c(100, 2), 
                options = list(x1D_stop = TRUE))
  
  # Try to run
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "afni_test_3ddecon")
  
  result <- tryCatch({
    run.afni_lm_spec(alm, outdir = test_dir, execute = TRUE)
    TRUE
  }, error = function(e) {
    message("3dDeconvolve error: ", e$message)
    FALSE
  })
  
  expect_true(result, "3dDeconvolve should execute without errors")
  
  # Check for output files
  if (result && dir.exists(test_dir)) {
    files <- list.files(test_dir)
    expect_gt(length(files), 0, "Should create output files")
    
    # Look for design matrix
    xmat_files <- files[grepl("xmat\\.1D$", files)]
    expect_gt(length(xmat_files), 0, "Should create design matrix file")
    
    # Clean up
    unlink(test_dir, recursive = TRUE)
  }
})

test_that("GLT commands are included when contrasts present", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # Should have GLT in command
  expect_equal(get_glt_count(cmd$cmd), 1)
  expect_true(grepl("-gltsym", cmd$cmd))
  expect_true(grepl("A_vs_B", cmd$cmd))
})

test_that("polort option controls polynomial detrending", {
  design <- make_simple_design(n_cond = 1, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  # Test different polort values
  cmd_neg1 <- generate_afni_command(model, options = list(polort = -1))
  expect_true(grepl("-polort -1", cmd_neg1$cmd))
  
  cmd_3 <- generate_afni_command(model, options = list(polort = 3))
  expect_true(grepl("-polort 3", cmd_3$cmd))
})

test_that("censor file and mask options work correctly", {
  design <- make_simple_design(n_cond = 2, n_runs = 2)
  sframe <- make_sampling_frame(n_runs = 2)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  # Create a mock censor vector (1=include, 0=exclude)
  censor <- c(rep(1, 90), rep(0, 10), rep(1, 100))  # Censor 10 TRs
  
  # Create dummy dataset with mask
  dset <- model$dataset
  dset$mask_file <- "brain_mask.nii"
  
  # Test with censor
  alm <- afni_lm(model, dset, 
    nodata = c(200, 2),
    censor = censor)
  
  # Verify censor file reference in command
  expect_true(grepl("-censor censor.1D", alm$cmd$cmd))
  
  # Check that censor is stored for file writing
  expect_equal(alm$cmd$censor, censor)
  
  # Test with mask
  alm_mask <- afni_lm(model, dset, nodata = c(200, 2))
  expect_true(grepl("-mask brain_mask.nii", alm_mask$cmd$cmd))
})

test_that("IRESP output and TENT basis function work together", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # TENT basis for estimating impulse response at multiple time points
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "TENT", start = 0, stop = 18, num = 10),
    sframe
  )
  
  # Generate command with IRESP output enabled
  dset <- model$dataset
  alm <- afni_lm(model, dset, 
    nodata = c(100, 2),
    options = list(iresp = TRUE))
  
  cmd <- alm$cmd
  
  # Check TENT parameters in command
  expect_true(grepl("TENT\\(0,18,10\\)", cmd$cmd))
  
  # Check IRESP output files are specified
  expect_true(grepl("-iresp", cmd$cmd))
  
  # Verify nbasis calculation for TENT (should be 10)
  term <- terms(model$event_model)[[1]]
  # TENT basis should have num basis functions
  expect_equal(attr(term$hrfspec, "num"), 10)
})

test_that("amplitude modulated design with stim_times_IM works", {
  # Create design with trial-wise amplitude modulation
  design <- data.frame(
    onset = seq(10, 190, by = 9),
    amplitude = rnorm(20, mean = 1, sd = 0.3),  # Trial-wise modulation
    run = 1
  )
  
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Use afni_trialwise for amplitude modulation
  emodel <- event_model(
    onset ~ afni_trialwise(amplitude, basis = "spmg1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Create dummy dataset and model
  dset <- list(scans = NULL, mask_file = NULL)
  bmodel <- baseline_model(basis = "constant", degree = 0, sframe = sframe)
  model <- fmri_model(emodel, bmodel, dset)
  
  alm <- afni_lm(model, dset, nodata = c(100, 2))
  cmd <- alm$cmd
  
  # Check for stim_times_IM instead of regular stim_times
  expect_true(grepl("-stim_times_IM", cmd$cmd))
  
  # Verify amplitude modulation is encoded in stimulus
  stim <- cmd$afni_stims[[1]]
  expect_s3_class(stim, "afni_stim_im_times")
  
  # Check that onsets include amplitude modulation (format: "onset*amplitude")
  # The onsets should be strings with asterisks for modulation
  expect_true(is.character(stim$onsets[[1]]))
  expect_true(grepl("\\*", stim$onsets[[1]]))
})