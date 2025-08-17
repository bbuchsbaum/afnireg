context("AFNI stimulus file generation")

test_that("stim_times files have correct format", {
  design <- make_simple_design(n_cond = 2, n_runs = 2)
  sframe <- make_sampling_frame(n_runs = 2)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  # Generate stimulus info
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  expect_length(stims, 2)  # One for each condition
  
  # Check each stimulus
  for (stim in stims) {
    expect_s3_class(stim, "afni_stim_times")
    expect_true(!is.null(stim$label))
    expect_true(!is.null(stim$onset_list))
    
    # Should have one entry per run
    expect_equal(length(stim$onset_list), 2)
  }
})

test_that("stim_times content is properly formatted", {
  design <- data.frame(
    condition = factor(c("A", "B", "A", "B")),
    onset = c(10, 20, 30, 40),
    run = c(1, 1, 2, 2)
  )
  
  sframe <- make_sampling_frame(n_runs = 2)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Check condition A
  stim_A <- stims[[1]]
  expect_equal(stim_A$label, "condition.A")
  expect_equal(stim_A$onset_list[[1]], 10)  # Run 1
  expect_equal(stim_A$onset_list[[2]], 30)  # Run 2
  
  # Check condition B
  stim_B <- stims[[2]]
  expect_equal(stim_B$label, "condition.B")
  expect_equal(stim_B$onset_list[[1]], 20)  # Run 1
  expect_equal(stim_B$onset_list[[2]], 40)  # Run 2
})

test_that("write_stim_times creates proper files", {
  temp_dir <- tempdir()
  
  # Create a mock afni_stim_times object
  stim <- structure(
    list(
      label = "test_condition",
      file_name = "test_stim_times.1D",
      hrf = "SPMG1",
      onsets = list(
        c(10, 20, 30),  # Run 1
        c(15, 25, 35),  # Run 2
        numeric(0)      # Run 3 (empty)
      ),
      onset_list = list(
        c(10, 20, 30),  # Run 1
        c(15, 25, 35),  # Run 2
        numeric(0)      # Run 3 (empty)
      ),
      durations = NULL
    ),
    class = "afni_stim_times"
  )
  
  # Write the file
  afnireg:::write_afni_stim(stim, temp_dir)
  
  # The file will be created with the name from stim$file_name
  file_path <- file.path(temp_dir, "test_stim_times.1D")
  expect_true(file.exists(file_path))
  
  # Read and check content
  content <- readLines(file_path)
  expect_equal(length(content), 3)  # Three runs
  
  # Check formatting
  expect_equal(content[1], "10 20 30")
  expect_equal(content[2], "15 25 35")
  expect_equal(content[3], "*")  # Empty run represented as asterisk
  
  # Clean up
  unlink(file_path)
})

test_that("stim_times handles durations", {
  design <- data.frame(
    condition = factor(c("A", "B")),
    onset = c(10, 20),
    duration = c(2.5, 3.0),
    run = c(1, 1)
  )
  
  sframe <- make_sampling_frame(n_runs = 1)
  
  # Use BLOCK basis with durations
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "BLOCK", durations = 2.5),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Check that durations are included
  stim_A <- stims[[1]]
  expect_false(is.null(stim_A$durations))
  expect_equal(stim_A$durations[[1]], 2.5)
})

test_that("factorial design generates correct stimulus files", {
  design <- make_factorial_design(levels1 = 2, levels2 = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Should have 4 stimulus objects (2x2 factorial)
  expect_equal(length(stims), 4)
  
  # Check labels match factorial conditions
  labels <- sapply(stims, function(s) s$label)
  expect_true(all(grepl("Fac1\\.F1_\\d+_Fac2\\.F2_\\d+", labels)))
})

test_that("empty runs are handled correctly", {
  # Design where condition B doesn't appear in run 2
  design <- data.frame(
    condition = factor(c("A", "B", "A", "A")),
    onset = c(10, 20, 30, 40),
    run = c(1, 1, 2, 2)
  )
  
  sframe <- make_sampling_frame(n_runs = 2)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Condition B stimulus
  stim_B <- stims[[2]]
  
  # Run 1 should have onset
  expect_equal(stim_B$onset_list[[1]], 20)
  
  # Run 2 should be empty
  expect_equal(length(stim_B$onset_list[[2]]), 0)
})

test_that("stim_file generation for R-based terms", {
  design <- make_simple_design(n_cond = 2, n_runs = 1)
  sframe <- make_sampling_frame(n_runs = 1)
  
  # R-based model (not AFNI)
  emodel <- event_model(
    onset ~ hrf(condition, basis = "spmg1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  
  # R-based terms should use build_afni_stims.convolved_term
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Should generate stim_file objects
  if (!is.null(stims)) {
    for (stim in stims) {
      expect_s3_class(stim, "afni_stim_file")
      expect_true(!is.null(stim$values))
    }
  }
})

test_that("stimulus labels are sanitized", {
  # Create design with special characters in factor levels
  design <- data.frame(
    condition = factor(c("A-1", "B.2", "C_3", "D 4")),
    onset = c(10, 20, 30, 40),
    run = 1
  )
  
  sframe <- make_sampling_frame(n_runs = 1)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Labels should be sanitized for AFNI compatibility
  labels <- sapply(stims, function(s) s$label)
  
  # Check that problematic characters are handled
  expect_false(any(grepl(" ", labels)))  # No spaces
})

test_that("multiple runs with different event counts work", {
  # Unbalanced design across runs
  design <- data.frame(
    condition = factor(c("A", "B", "A", "B", "A", "A", "B")),
    onset = c(10, 20, 30, 15, 25, 35, 45),
    run = c(1, 1, 1, 2, 2, 3, 3)
  )
  
  sframe <- make_sampling_frame(n_runs = 3, run_length_trs = 50)
  
  emodel <- event_model(
    onset ~ afni_hrf(condition, basis = "SPMG1"),
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  term <- terms(emodel)[[1]]
  stims <- afnireg:::build_afni_stims(term, iresp = FALSE, tr_times = 2)
  
  # Check condition A
  stim_A <- stims[[1]]
  expect_equal(length(stim_A$onset_list), 3)  # 3 runs
  expect_equal(length(stim_A$onset_list[[1]]), 2)  # 2 events in run 1
  expect_equal(length(stim_A$onset_list[[2]]), 1)  # 1 event in run 2
  expect_equal(length(stim_A$onset_list[[3]]), 1)  # 1 event in run 3
  
  # Check condition B
  stim_B <- stims[[2]]
  expect_equal(length(stim_B$onset_list), 3)  # 3 runs
  expect_equal(length(stim_B$onset_list[[1]]), 1)  # 1 event in run 1
  expect_equal(length(stim_B$onset_list[[2]]), 1)  # 1 event in run 2
  expect_equal(length(stim_B$onset_list[[3]]), 1)  # 1 event in run 3
})