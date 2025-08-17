context("AFNI mixed models (R and AFNI terms)")

test_that("mixed model with hrf() and afni_hrf() creates proper terms", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Create mixed model with both R-based and AFNI-native terms
  emodel <- event_model(
    onset ~ hrf(RT, basis = "spmg1") + afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Check we have two terms
  expect_length(terms(emodel), 2)
  
  # First term should be R-based (event_term)
  expect_s3_class(terms(emodel)[[1]], "event_term")
  expect_false(requires_external_processing(terms(emodel)[[1]]))
  
  # Second term should be AFNI-based
  expect_s3_class(terms(emodel)[[2]], "afni_hrf_convolved_term")
  expect_true(requires_external_processing(terms(emodel)[[2]]))
})

test_that("R-based terms generate stim_file commands", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Pure R model
  model <- build_test_afni_model(
    design,
    onset ~ hrf(condition, basis = "spmg1"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # R-based terms should use -stim_file
  expect_has_stim_file(cmd$cmd, 2)  # One for each condition
  expect_has_stim_times(cmd$cmd, 0)  # No -stim_times
})

test_that("AFNI terms generate stim_times commands", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Pure AFNI model
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # AFNI terms should use -stim_times
  expect_has_stim_times(cmd$cmd, 2)  # One for each condition
  expect_has_stim_file(cmd$cmd, 0)  # No -stim_file
})

test_that("mixed model generates both stim_file and stim_times", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Mixed model
  model <- build_test_afni_model(
    design,
    onset ~ hrf(RT, basis = "spmg1") + afni_hrf(condition, basis = "SPMG1"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # Should have both types of stims
  expect_has_stim_file(cmd$cmd, 1)   # One for RT (continuous)
  expect_has_stim_times(cmd$cmd, 2)  # Two for conditions A and B
})

test_that("external processing detection works correctly", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ hrf(RT) + afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Check each term
  r_term <- terms(emodel)[[1]]
  afni_term <- terms(emodel)[[2]]
  
  expect_false(requires_external_processing(r_term))
  expect_true(requires_external_processing(afni_term))
})

test_that("contrast_weights.event_model skips AFNI terms", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  # Mixed model with contrasts only on AFNI term
  emodel <- event_model(
    onset ~ hrf(run) + afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Try to get contrast weights from event_model
  # This should skip the AFNI term and not error
  weights <- tryCatch({
    contrast_weights(emodel)
  }, error = function(e) NULL)
  
  # Should not error
  expect_false(is.null(weights))
})

test_that("design matrix excludes AFNI terms", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Mixed model
  emodel <- event_model(
    onset ~ hrf(RT, basis = "spmg1") + afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  dm <- design_matrix(emodel)
  
  # Design matrix should only include R-based term columns
  # RT is continuous, so should have 1 column
  expect_equal(ncol(dm), 1)
  expect_true(grepl("RT", names(dm)[1]))
})

test_that("col_indices handles mixed models correctly", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ hrf(RT) + afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  dm <- design_matrix(emodel)
  col_idx <- attr(dm, "col_indices")
  
  expect_type(col_idx, "list")
  expect_length(col_idx, 2)
  
  # First term (R-based) should have indices
  expect_gt(length(col_idx[[1]]), 0)
  
  # Second term (AFNI) should have no indices
  expect_equal(length(col_idx[[2]]), 0)
})

test_that("contrasts work in mixed models", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  model <- build_test_afni_model(
    design,
    onset ~ hrf(RT) + afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # Should have GLT from AFNI term
  expect_gt(length(cmd$glts), 0)
  expect_true("A_vs_B" %in% cmd$gltnames)
})

test_that("multiple AFNI terms work", {
  design <- data.frame(
    cond1 = factor(rep(c("A", "B"), 10)),
    cond2 = factor(rep(c("X", "Y"), each = 10)),
    onset = seq(1, 200, length.out = 20),
    run = 1
  )
  
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Multiple AFNI terms
  emodel <- event_model(
    onset ~ afni_hrf(cond1, basis = "SPMG1") + afni_hrf(cond2, basis = "GAM"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Both should be AFNI terms
  expect_true(all(sapply(terms(emodel), requires_external_processing)))
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(cond1, basis = "SPMG1") + afni_hrf(cond2, basis = "GAM"),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # Should have stim_times for all conditions
  expect_has_stim_times(cmd$cmd, 4)  # A, B, X, Y
})

test_that("R and AFNI terms can have different basis functions", {
  design <- make_mixed_model_design(n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Different bases for R and AFNI
  emodel <- event_model(
    onset ~ hrf(RT, basis = "spmg3") + afni_hrf(condition, basis = "BLOCK"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # R term should have 3 basis functions (spmg3)
  # Access HRF directly from hrfspec attribute for event_term
  term1 <- terms(emodel)[[1]]
  hrfspec1 <- attr(term1, "hrfspec")
  expect_equal(attr(hrfspec1$hrf, "nbasis"), 3)
  
  # AFNI term should have 1 basis function
  expect_equal(nbasis(terms(emodel)[[2]]), 1)
})

test_that("AFNI registration persists across package loads", {
  # Check that AFNI terms are registered
  expect_true(is_external_hrfspec("afni_hrfspec"))
  expect_true(is_external_hrfspec("afni_hrf_convolved_term"))
  
  # Check registration info
  info <- get_external_hrfspec_info("afni_hrfspec")
  expect_equal(info$package, "afnireg")
  expect_true(info$requires_external_processing)
})