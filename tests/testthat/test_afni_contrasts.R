context("AFNI contrast handling")

test_that("afni_hrf preserves single contrast", {
  con <- contrast(~ A - B, name = "A_vs_B")
  spec <- afni_hrf(condition, basis = "SPMG1", contrasts = con)
  
  expect_length(spec$contrasts, 1)
  expect_equal(spec$contrasts[[1]]$name, "A_vs_B")
})

test_that("afni_hrf preserves contrast list", {
  con1 <- contrast(~ A - B, name = "A_vs_B")
  con2 <- contrast(~ A + B, name = "A_and_B")
  
  spec <- afni_hrf(condition, basis = "SPMG1", contrasts = list(con1, con2))
  
  expect_length(spec$contrasts, 2)
  expect_equal(spec$contrasts[[1]]$name, "A_vs_B")
  expect_equal(spec$contrasts[[2]]$name, "A_and_B")
})

test_that("afni_hrf preserves contrast_set", {
  con1 <- contrast(~ A - B, name = "A_vs_B")
  con2 <- contrast(~ A + B, name = "A_and_B")
  cset <- contrast_set(con1, con2)
  
  spec <- afni_hrf(condition, basis = "SPMG1", contrasts = cset)
  
  expect_length(spec$contrasts, 2)
  expect_s3_class(spec$contrasts, "contrast_set")
})

test_that("contrast_weights.afni_hrf_convolved_term returns weights", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  weights <- contrast_weights(term)
  
  expect_type(weights, "list")
  expect_length(weights, 1)
  expect_equal(names(weights), "A_vs_B")
  
  # Check weight values
  w <- weights[[1]]$weights
  expect_equal(as.vector(w), c(1, -1))
})

test_that("contrast_weights handles factorial design contrasts", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Main effect contrast
  con_main <- contrast(
    ~ Fac1.F1_1_Fac2.F2_1 + Fac1.F1_1_Fac2.F2_2 - 
      Fac1.F1_2_Fac2.F2_1 - Fac1.F1_2_Fac2.F2_2,
    name = "Main_Fac1"
  )
  
  emodel <- event_model(
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = list(con_main)),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  weights <- contrast_weights(term)
  
  expect_length(weights, 1)
  expect_equal(names(weights), "Main_Fac1")
  
  w <- weights[[1]]$weights
  # Main effect of Fac1: F1_1 levels get +1, F1_2 levels get -1
  # Order follows longnames: F1_1/F2_1, F1_2/F2_1, F1_1/F2_2, F1_2/F2_2
  expect_equal(as.vector(w), c(1, -1, 1, -1))
})

test_that("to_glt converts contrast weights to GLT format", {
  # Create a simple contrast object
  con_obj <- list(
    weights = c(1, -1, 0, 0),
    condnames = c("A", "B", "C", "D"),
    name = "A_vs_B"
  )
  class(con_obj) <- "contrast"
  
  glt <- to_glt(con_obj)
  
  expect_s3_class(glt, "glt")
  expect_equal(glt$name, "A_vs_B")
  expect_equal(glt$glt_str, "1*A -1*B")
})

test_that("to_glt handles complex contrasts", {
  con_obj <- list(
    weights = c(0.5, 0.5, -0.5, -0.5),
    condnames = c("A1", "A2", "B1", "B2"),
    name = "A_vs_B"
  )
  class(con_obj) <- "contrast"
  
  glt <- to_glt(con_obj)
  
  expect_equal(glt$glt_str, "0.5*A1 +0.5*A2 -0.5*B1 -0.5*B2")
})

test_that("build_decon_command extracts contrasts from AFNI terms", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  # Check that contrasts were extracted
  expect_gt(length(cmd$glts), 0)
  expect_true("A_vs_B" %in% cmd$gltnames)
})

test_that("multiple contrasts generate multiple GLTs", {
  design <- make_simple_design(n_cond = 3, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con1 <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  con2 <- contrast(~ condition.A - condition.C, name = "A_vs_C")
  con3 <- contrast(~ condition.B - condition.C, name = "B_vs_C")
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(condition, basis = "SPMG1", 
                    contrasts = list(con1, con2, con3)),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  expect_equal(length(cmd$glts), 3)
  expect_equal(sort(cmd$gltnames), sort(c("A_vs_B", "A_vs_C", "B_vs_C")))
})

test_that("factorial contrasts work correctly", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Create all standard contrasts for 2x2 factorial
  con_main1 <- contrast(
    ~ Fac1.F1_1_Fac2.F2_1 + Fac1.F1_1_Fac2.F2_2 - 
      Fac1.F1_2_Fac2.F2_1 - Fac1.F1_2_Fac2.F2_2,
    name = "Main_Fac1"
  )
  
  con_main2 <- contrast(
    ~ Fac1.F1_1_Fac2.F2_1 + Fac1.F1_2_Fac2.F2_1 - 
      Fac1.F1_1_Fac2.F2_2 - Fac1.F1_2_Fac2.F2_2,
    name = "Main_Fac2"
  )
  
  con_int <- contrast(
    ~ (Fac1.F1_1_Fac2.F2_1 - Fac1.F1_1_Fac2.F2_2) - 
      (Fac1.F1_2_Fac2.F2_1 - Fac1.F1_2_Fac2.F2_2),
    name = "Interaction"
  )
  
  all_contrasts <- contrast_set(con_main1, con_main2, con_int)
  
  model <- build_test_afni_model(
    design,
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts),
    sframe
  )
  
  cmd <- generate_afni_command(model)
  
  expect_equal(length(cmd$glts), 3)
  expect_true("Main_Fac1" %in% cmd$gltnames)
  expect_true("Main_Fac2" %in% cmd$gltnames)
  expect_true("Interaction" %in% cmd$gltnames)
})

test_that("contrast extraction works with construct.afni_hrfspec", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ A - B, name = "A_vs_B")
  spec <- afni_hrf(condition, basis = "SPMG1", contrasts = list(con))
  
  model_spec <- list(
    onsets = design$onset,
    blockids = design$run,
    durations = rep(0, nrow(design)),
    data = design,
    formula_env = environment(),
    sampling_frame = sframe
  )
  
  term <- construct(spec, model_spec)
  
  # Check that contrasts are preserved in the constructed term
  hrfspec_after <- attr(term, "hrfspec")
  expect_equal(length(hrfspec_after$contrasts), 1)
  expect_equal(hrfspec_after$contrasts[[1]]$name, "A_vs_B")
})

test_that("contrasts method returns preserved contrasts", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  con <- contrast(~ condition.A - condition.B, name = "A_vs_B")
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1", contrasts = list(con)),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  term_contrasts <- contrasts(term)
  
  expect_length(term_contrasts, 1)
  expect_equal(term_contrasts[[1]]$name, "A_vs_B")
})

test_that("GLT string formatting is correct", {
  # Test positive weights
  con_obj <- list(
    weights = c(1, 0, 0),
    condnames = c("A", "B", "C"),
    name = "A_only"
  )
  class(con_obj) <- "contrast"
  
  glt <- to_glt(con_obj)
  expect_equal(glt$glt_str, "1*A")
  
  # Test negative weights
  con_obj$weights <- c(-1, 0, 0)
  glt <- to_glt(con_obj)
  expect_equal(glt$glt_str, "-1*A")
  
  # Test fractional weights
  con_obj$weights <- c(0.5, -0.25, 0.25)
  glt <- to_glt(con_obj)
  expect_equal(glt$glt_str, "0.5*A -0.25*B +0.25*C")
})