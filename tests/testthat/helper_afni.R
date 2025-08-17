# Helper functions for AFNI tests

library(fmridesign)
library(fmrihrf)
library(afnireg)

# Check if 3dDeconvolve is available
has_3ddeconvolve <- function() {
  file.exists("~/abin/3dDeconvolve") || 
  file.exists("/usr/local/abin/3dDeconvolve") ||
  Sys.which("3dDeconvolve") != ""
}

# Create simple event design
make_simple_design <- function(n_cond = 2, n_runs = 1, n_events_per_cond = 4) {
  if (n_cond == 1) {
    # Single condition
    design <- data.frame(
      condition = factor(rep("A", n_events_per_cond * n_runs)),
      run = rep(1:n_runs, each = n_events_per_cond)
    )
  } else {
    # Multiple conditions
    conditions <- LETTERS[1:n_cond]
    design <- expand.grid(
      condition = factor(conditions),
      rep = 1:n_events_per_cond,
      run = 1:n_runs
    )
  }
  
  # Generate onsets
  run_length <- 100  # seconds
  for (r in 1:n_runs) {
    run_idx <- which(design$run == r)
    n_events <- length(run_idx)
    design$onset[run_idx] <- seq(10, run_length - 10, length.out = n_events)
  }
  
  design
}

# Create factorial design
make_factorial_design <- function(levels1 = 2, levels2 = 2, n_reps = 4, n_runs = 1) {
  fac1_levels <- paste0("F1_", 1:levels1)
  fac2_levels <- paste0("F2_", 1:levels2)
  
  design <- expand.grid(
    Fac1 = factor(fac1_levels),
    Fac2 = factor(fac2_levels),
    rep = 1:n_reps,
    run = 1:n_runs
  )
  
  # Generate onsets
  run_length <- 200  # seconds
  for (r in 1:n_runs) {
    run_idx <- which(design$run == r)
    n_events <- length(run_idx)
    base_onsets <- seq(10, run_length - 10, length.out = n_events)
    # Add small jitter
    design$onset[run_idx] <- sort(base_onsets + runif(n_events, -2, 2))
  }
  
  design
}

# Create mixed model design with categorical and continuous variables
make_mixed_model_design <- function(n_runs = 1) {
  design <- data.frame(
    condition = factor(rep(c("A", "B"), each = 10)),
    RT = rnorm(20, mean = 0.5, sd = 0.2),
    run = rep(1:n_runs, length.out = 20)
  )
  
  design$onset <- seq(1, 200, length.out = 20)
  design
}

# Create standard sampling frame
make_sampling_frame <- function(n_runs = 1, run_length_trs = 100, TR = 2) {
  sampling_frame(blocklens = rep(run_length_trs, n_runs), TR = TR)
}

# Create dummy dataset for testing
make_dummy_dataset <- function(event_table, TR = 2, run_length = 100) {
  fmri_dataset(
    scans = "dummy.nii",
    mask = "dummy_mask.nii",
    TR = TR,
    run_length = run_length,
    event_table = event_table,
    base_path = ".",
    dummy_mode = TRUE
  )
}

# Check if command contains expected stim_times entries
expect_has_stim_times <- function(cmd, n_expected, info = NULL) {
  matches <- gregexpr("-stim_times", cmd)[[1]]
  n_found <- if (matches[1] == -1) 0 else length(matches)
  expect_equal(n_found, n_expected, info = info)
}

# Check if command contains expected stim_file entries
expect_has_stim_file <- function(cmd, n_expected, info = NULL) {
  matches <- gregexpr("-stim_file", cmd)[[1]]
  n_found <- if (matches[1] == -1) 0 else length(matches)
  expect_equal(n_found, n_expected, info = info)
}

# Check if command contains expected GLT entries
expect_has_glt <- function(cmd, glt_names, info = NULL) {
  for (name in glt_names) {
    expect_true(grepl(name, cmd), 
                info = paste0("Expected GLT '", name, "' not found. ", info))
  }
}

# Check if command contains nodata flag
expect_has_nodata <- function(cmd, NT, TR, info = NULL) {
  pattern <- sprintf("-nodata %d %d", NT, TR)
  expect_true(grepl(pattern, cmd), 
              info = paste0("Expected '", pattern, "' not found. ", info))
}

# Extract GLT count from command
get_glt_count <- function(cmd) {
  match <- regexpr("-num_glt (\\d+)", cmd)
  if (match > 0) {
    matched <- regmatches(cmd, match)
    as.integer(sub("-num_glt ", "", matched))
  } else {
    0
  }
}

# Extract stim count from command
get_stim_count <- function(cmd) {
  match <- regexpr("-num_stimts (\\d+)", cmd)
  if (match > 0) {
    matched <- regmatches(cmd, match)
    as.integer(sub("-num_stimts ", "", matched))
  } else {
    0
  }
}

# Build a standard AFNI model for testing
build_test_afni_model <- function(design, formula, sframe = NULL, 
                                 baseline_degree = 2, contrasts = NULL) {
  if (is.null(sframe)) {
    n_runs <- length(unique(design$run))
    sframe <- make_sampling_frame(n_runs)
  }
  
  # Create event model
  emodel <- event_model(
    formula,
    data = design,
    block = ~ run,
    sampling_frame = sframe
  )
  
  # Create baseline model
  bmodel <- baseline_model(
    basis = if (baseline_degree == 0) "constant" else "poly",
    degree = baseline_degree,
    sframe = sframe
  )
  
  # Create dummy dataset
  dset <- make_dummy_dataset(design)
  
  # Create full model
  fmri_model(emodel, bmodel, dset)
}

# Generate AFNI command from model
generate_afni_command <- function(model, nodata = c(100, 2), options = list()) {
  dset <- model$dataset
  if (is.null(dset)) {
    dset <- list(scans = NULL, mask_file = NULL)
  }
  
  default_opts <- list(
    nodata = nodata,
    polort = -1,
    iresp = FALSE,
    TR_times = nodata[2],
    x1D_stop = TRUE,
    nofullf_atall = TRUE
  )
  
  opts <- modifyList(default_opts, options)
  
  afnireg:::build_decon_command(
    model = model,
    dataset = dset,
    working_dir = tempdir(),
    opts = opts
  )
}

# Check stimulus file content (for debugging)
check_stim_file_content <- function(file_path, expected_runs = NULL) {
  if (!file.exists(file_path)) {
    return(list(exists = FALSE))
  }
  
  content <- readLines(file_path)
  
  result <- list(
    exists = TRUE,
    n_lines = length(content),
    content = content
  )
  
  if (!is.null(expected_runs)) {
    expect_equal(length(content), expected_runs,
                info = sprintf("File %s should have %d lines (one per run)", 
                              file_path, expected_runs))
  }
  
  result
}