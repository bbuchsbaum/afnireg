context("AFNI HRF core functionality")

test_that("afni_hrf creates proper hrfspec object", {
  # Test basic creation
  spec <- afni_hrf(condition, basis = "SPMG1")
  
  expect_s3_class(spec, "afni_hrfspec")
  expect_s3_class(spec, "hrfspec")
  expect_equal(as.vector(spec$hrf), "SPMG1")
  expect_equal(spec$name, "condition")  # name is the variable name, not "afni_hrf"
})

test_that("afni_hrf accepts various AFNI basis functions", {
  # Only test bases that are actually supported (check match.arg in afni_hrf)
  bases <- c("spmg1", "gam", "block", "sin", "tent")
  expected <- c("SPMG1", "GAM", "BLOCK", "SIN", "TENT")
  
  for (i in seq_along(bases)) {
    spec <- afni_hrf(condition, basis = bases[i])
    expect_equal(as.vector(spec$hrf), expected[i],
                info = paste("Failed for basis:", bases[i]))
  }
})

test_that("afni_hrf accepts lowercase aliases", {
  spec_sin <- afni_hrf(condition, basis = "sin")
  spec_gam <- afni_hrf(condition, basis = "gam")
  spec_block <- afni_hrf(condition, basis = "block")
  
  expect_equal(as.vector(spec_sin$hrf), "SIN")
  expect_equal(as.vector(spec_gam$hrf), "GAM")
  expect_equal(as.vector(spec_block$hrf), "BLOCK")
})

test_that("afni_hrf preserves contrasts", {
  con1 <- contrast(~ A - B, name = "A_vs_B")
  con2 <- contrast(~ A + B, name = "A_and_B")
  contrast_list <- list(con1, con2)
  
  spec <- afni_hrf(condition, basis = "SPMG1", contrasts = contrast_list)
  
  expect_equal(length(spec$contrasts), 2)
  expect_equal(spec$contrasts[[1]]$name, "A_vs_B")
  expect_equal(spec$contrasts[[2]]$name, "A_and_B")
})

test_that("construct.afni_hrfspec creates afni_hrf_convolved_term", {
  # Create a simple design
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Create hrfspec
  spec <- afni_hrf(condition, basis = "SPMG1")
  
  # Create model spec for construct
  model_spec <- list(
    onsets = design$onset,
    blockids = design$run,
    durations = rep(0, nrow(design)),
    data = design,
    formula_env = environment(),
    sampling_frame = sframe
  )
  
  # Construct the term
  term <- construct(spec, model_spec)
  
  expect_s3_class(term, "afni_hrf_convolved_term")
  expect_s3_class(term, "event_term")
  expect_equal(attr(term, "hrfspec"), spec)
})

test_that("longnames returns correct names for simple design", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Create model with afni_hrf
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  lnames <- longnames(term)
  
  expect_equal(length(lnames), 2)
  expect_true("condition.A" %in% lnames)
  expect_true("condition.B" %in% lnames)
})

test_that("longnames handles factorial designs", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Create model with factorial afni_hrf
  emodel <- event_model(
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  lnames <- longnames(term)
  
  expect_equal(length(lnames), 4)  # 2x2 factorial = 4 conditions
  expect_true(all(grepl("Fac1\\.F1_\\d+_Fac2\\.F2_\\d+", lnames)))
})

test_that("shortnames matches longnames for factorial designs", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  # For AFNI terms with factorial designs, shortnames should match longnames
  expect_equal(shortnames(term), longnames(term))
})

test_that("nbasis returns correct value for AFNI HRF", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  # AFNI HRFs have single basis function
  expect_equal(nbasis(term), 1)
})

test_that("is_continuous returns FALSE for AFNI categorical terms", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  expect_false(is_continuous(term))
})

test_that("cells returns correct factor combinations", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  cell_df <- cells(term)
  
  expect_equal(nrow(cell_df), 4)  # 2x2 factorial
  expect_equal(ncol(cell_df), 2)  # 2 factors
  expect_equal(names(cell_df), c("Fac1", "Fac2"))
})

test_that("conditions returns proper condition names", {
  design <- make_simple_design(n_cond = 3, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  conds <- conditions(term)
  
  expect_equal(length(conds), 3)
  expect_true(all(c("condition.A", "condition.B", "condition.C") %in% conds))
})

test_that("split_onsets handles simple designs correctly", {
  design <- make_simple_design(n_cond = 2, n_runs = 2)
  sframe <- make_sampling_frame(n_runs = 2)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  split_result <- fmridesign::split_onsets(term$evterm, sframe, global=FALSE, blocksplit = TRUE)
  
  expect_type(split_result, "list")
  expect_equal(length(split_result), 2)  # 2 conditions
  
  # Each condition should have onsets for each run
  for (onsets_list in split_result) {
    expect_equal(length(onsets_list), 2)  # 2 runs
  }
})

test_that("AFNI terms are registered as requiring external processing", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  expect_true(requires_external_processing(term))
})

test_that("afni_hrf handles durations parameter", {
  spec <- afni_hrf(condition, basis = "BLOCK", durations = 2.5)
  
  expect_equal(spec$durations, 2.5)
})

test_that("afni_hrf handles subset parameter", {
  spec <- afni_hrf(condition, basis = "SPMG1", subset = condition == "A")
  
  expect_false(is.null(spec$subset))
})