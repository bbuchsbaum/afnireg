#' AFNI Test Suite
#' 
#' Systematic tests for AFNI integration, building from simple to complex

library(fmrireg)
library(fmridesign)
library(fmrihrf)
library(afnireg)

# Source the test runner
source("/Users/bbuchsbaum/code/afnireg/afni_test/run_afni_test.R")

# Store all test results
all_tests <- list()

# ==============================================================================
# TEST 1: Simple single stimulus, single run
# ==============================================================================
run_test_1 <- function() {
  message("\n" , paste(rep("=", 70), collapse=""))
  message("TEST 1: Simple Single Stimulus")
  message(paste(rep("=", 70), collapse=""))
  
  # Create minimal event table with one stimulus
  design <- data.frame(
    stim = factor("A"),
    onset = 10,  # Single event at 10 seconds
    run = 1
  )
  
  # Set up sampling frame (100 time points, TR=2)
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Minimal baseline (just intercept/constant)
  bmodel <- baseline_model(basis = "constant", sframe = sframe)
  
  # Run test
  result <- run_afni_test(
    test_name = "test1_simple_single_stim",
    design = design,
    model_formula = onset ~ afni_hrf(stim, basis = "SPMG1"),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "Simplest case: one stimulus type, one event, one run",
    nodata = c(100, 2)
  )
  
  return(result)
}

# ==============================================================================
# TEST 2: Two conditions (A and B), single run
# ==============================================================================
run_test_2 <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST 2: Two Conditions Single Run")
  message(paste(rep("=", 70), collapse=""))
  
  # Create event table with two conditions
  design <- data.frame(
    stim = factor(c("A", "B", "A", "B")),
    onset = c(10, 30, 50, 70),
    run = 1
  )
  
  # Set up sampling frame
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Baseline with polynomial drift correction
  bmodel <- baseline_model(basis = "poly", degree = 2, sframe = sframe)
  
  # Run test
  result <- run_afni_test(
    test_name = "test2_two_conditions",
    design = design,
    model_formula = onset ~ afni_hrf(stim, basis = "SPMG1"),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "Two stimulus conditions (A/B) with polynomial baseline",
    nodata = c(100, 2),
    options = list(x1D_stop = TRUE)
  )
  
  return(result)
}

# ==============================================================================
# TEST 3: Multiple runs
# ==============================================================================
run_test_3 <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST 3: Multiple Runs")
  message(paste(rep("=", 70), collapse=""))
  
  # Create event table with multiple runs
  design <- data.frame(
    stim = factor(c("A", "B", "A", "B", "A", "B")),
    onset = c(10, 30, 50, 10, 30, 50),
    run = c(1, 1, 1, 2, 2, 2)
  )
  
  # Set up sampling frame with 2 runs
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  
  # Baseline model
  bmodel <- baseline_model(basis = "poly", degree = 1, sframe = sframe)
  
  # Run test
  result <- run_afni_test(
    test_name = "test3_multiple_runs",
    design = design,
    model_formula = onset ~ afni_hrf(stim, basis = "GAM"),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "Two runs with A/B conditions, testing run concatenation",
    nodata = c(200, 2)
  )
  
  return(result)
}

# ==============================================================================
# TEST 4: Different AFNI HRF models
# ==============================================================================
run_test_4 <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST 4: Different HRF Models")
  message(paste(rep("=", 70), collapse=""))
  
  # Create event table
  design <- data.frame(
    cond = factor(c("task", "control", "task", "control")),
    onset = c(10, 40, 70, 100),
    run = 1
  )
  
  # Set up sampling frame
  sframe <- sampling_frame(blocklens = 150, TR = 2)
  
  # Baseline model
  bmodel <- baseline_model(basis = "constant", sframe = sframe)
  
  # Test with BLOCK HRF
  result <- run_afni_test(
    test_name = "test4_block_hrf",
    design = design,
    model_formula = onset ~ afni_hrf(cond, basis = "BLOCK(10,1)"),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "Testing BLOCK HRF model for block design",
    nodata = c(150, 2)
  )
  
  return(result)
}

# ==============================================================================
# TEST 5: With duration information
# ==============================================================================
run_test_5 <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST 5: Events with Durations")
  message(paste(rep("=", 70), collapse=""))
  
  # Create event table with durations
  design <- data.frame(
    stim = factor(c("short", "long", "short", "long")),
    onset = c(10, 30, 60, 90),
    duration = c(2, 10, 2, 10),
    run = 1
  )
  
  # Set up sampling frame
  sframe <- sampling_frame(blocklens = 150, TR = 2)
  
  # Baseline model
  bmodel <- baseline_model(basis = "poly", degree = 2, sframe = sframe)
  
  # Run test with duration-modulated HRF
  result <- run_afni_test(
    test_name = "test5_durations",
    design = design,
    model_formula = onset ~ afni_hrf(stim, basis = "dmBLOCK"),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "Events with varying durations using dmBLOCK",
    nodata = c(150, 2)
  )
  
  return(result)
}

# ==============================================================================
# TEST 6: 2x2 Factorial Design with Contrasts
# ==============================================================================
run_test_6 <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST 6: 2x2 Factorial Design with Contrasts")
  message(paste(rep("=", 70), collapse=""))
  
  # Create 2x2 factorial design with 3 runs
  # Each condition appears 4 times per run
  # This creates all 4 combinations: aa-cc, aa-dd, bb-cc, bb-dd
  design <- expand.grid(
    Fac1 = factor(c("aa", "bb")),
    Fac2 = factor(c("cc", "dd")),
    rep = 1:4,
    run = 1:3
  )
  
  # Generate onsets (spread events across each run)
  # Each run is 100 TRs (200 seconds) with TR=2
  run_length_trs <- 100  # Length of each run in TRs
  run_length_secs <- run_length_trs * 2  # Length in seconds
  events_per_run <- 16  # 4 conditions x 4 repetitions
  
  design <- design[order(design$run, design$rep, design$Fac1, design$Fac2), ]
  design$onset <- NA
  
  for (r in 1:3) {
    run_idx <- which(design$run == r)
    # Spread events evenly across the run with some jitter
    base_onsets <- seq(10, run_length_secs - 10, length.out = events_per_run)
    # Add small jitter to avoid perfect periodicity
    jittered_onsets <- base_onsets + runif(events_per_run, -2, 2)
    design$onset[run_idx] <- sort(jittered_onsets)
  }
  
  # Set up sampling frame with 3 runs (blocklens in TRs)
  sframe <- sampling_frame(blocklens = c(run_length_trs, run_length_trs, run_length_trs), TR = 2)
  
  # Baseline model with polynomial drift
  bmodel <- baseline_model(basis = "poly", degree = 2, sframe = sframe)
  
  # Define contrasts using the actual condition names from the factorial design
  library(fmridesign)
  
  # Main effect of Fac1 (aa vs bb)
  con_fac1 <- contrast(
    ~ Fac1.aa_Fac2.cc + Fac1.aa_Fac2.dd - Fac1.bb_Fac2.cc - Fac1.bb_Fac2.dd,
    name = "Main_Fac1_aa_vs_bb"
  )
  
  # Main effect of Fac2 (cc vs dd)
  con_fac2 <- contrast(
    ~ Fac1.aa_Fac2.cc + Fac1.bb_Fac2.cc - Fac1.aa_Fac2.dd - Fac1.bb_Fac2.dd,
    name = "Main_Fac2_cc_vs_dd"
  )
  
  # Interaction contrast: (aa.cc - aa.dd) - (bb.cc - bb.dd)
  # This tests whether the effect of Fac2 differs between levels of Fac1
  con_interaction <- contrast(
    ~ (Fac1.aa_Fac2.cc - Fac1.aa_Fac2.dd) - (Fac1.bb_Fac2.cc - Fac1.bb_Fac2.dd),
    name = "Interaction_Fac1xFac2"
  )
  
  # All pairwise comparisons between cells using the actual condition names
  con_aa_cc_vs_aa_dd <- contrast(
    ~ Fac1.aa_Fac2.cc - Fac1.aa_Fac2.dd,
    name = "aa.cc_vs_aa.dd"
  )
  
  con_aa_cc_vs_bb_cc <- contrast(
    ~ Fac1.aa_Fac2.cc - Fac1.bb_Fac2.cc,
    name = "aa.cc_vs_bb.cc"
  )
  
  con_aa_cc_vs_bb_dd <- contrast(
    ~ Fac1.aa_Fac2.cc - Fac1.bb_Fac2.dd,
    name = "aa.cc_vs_bb.dd"
  )
  
  con_aa_dd_vs_bb_cc <- contrast(
    ~ Fac1.aa_Fac2.dd - Fac1.bb_Fac2.cc,
    name = "aa.dd_vs_bb.cc"
  )
  
  con_aa_dd_vs_bb_dd <- contrast(
    ~ Fac1.aa_Fac2.dd - Fac1.bb_Fac2.dd,
    name = "aa.dd_vs_bb.dd"
  )
  
  con_bb_cc_vs_bb_dd <- contrast(
    ~ Fac1.bb_Fac2.cc - Fac1.bb_Fac2.dd,
    name = "bb.cc_vs_bb.dd"
  )
  
  # Combine all contrasts using contrast_set
  all_contrasts <- do.call(contrast_set, list(
    con_fac1,
    con_fac2,
    con_interaction,
    con_aa_cc_vs_aa_dd,
    con_aa_cc_vs_bb_cc,
    con_aa_cc_vs_bb_dd,
    con_aa_dd_vs_bb_cc,
    con_aa_dd_vs_bb_dd,
    con_bb_cc_vs_bb_dd
  ))
  
  # Run test with contrasts now that the evaluation works
  result <- run_afni_test(
    test_name = "test6_factorial_2x2",
    design = design,
    model_formula = onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts),
    sampling_frame = sframe,
    baseline_model = bmodel,
    test_description = "2x2 factorial design (Fac1: aa/bb, Fac2: cc/dd) with 3 runs and all pairwise contrasts",
    nodata = c(run_length_trs * 3, 2),  # Total TRs across all 3 runs
    contrasts = NULL,
    options = list(x1D_stop = TRUE)
  )
  
  # Verify factorial design created correct stimulus files
  if (!is.null(result$test_spec$files_created)) {
    stim_files <- grep("_times\\.1D$", result$test_spec$files_created, value = TRUE)
    if (length(stim_files) == 4) {
      message("✓ Successfully created 4 stimulus files for 2x2 factorial design:")
      for (f in stim_files) {
        message("  - ", f)
      }
    } else {
      message("⚠ Warning: Expected 4 stimulus files, found ", length(stim_files))
    }
  }
  
  # Note: GLT conversion from contrasts requires additional implementation
  # The contrasts are properly defined but need formula evaluation context handling
  
  return(result)
}

# ==============================================================================
# RUN ALL TESTS
# ==============================================================================
run_all_tests <- function() {
  message("\n", paste(rep("=", 70), collapse=""))
  message("RUNNING AFNI TEST SUITE")
  message(paste(rep("=", 70), collapse=""))
  
  test_results <- list()
  
  # Run each test
  tryCatch({
    test_results$test1 <- run_test_1()
  }, error = function(e) {
    message("TEST 1 FAILED: ", e$message)
    test_results$test1 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  tryCatch({
    test_results$test2 <- run_test_2()
  }, error = function(e) {
    message("TEST 2 FAILED: ", e$message)
    test_results$test2 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  tryCatch({
    test_results$test3 <- run_test_3()
  }, error = function(e) {
    message("TEST 3 FAILED: ", e$message)
    test_results$test3 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  tryCatch({
    test_results$test4 <- run_test_4()
  }, error = function(e) {
    message("TEST 4 FAILED: ", e$message)
    test_results$test4 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  tryCatch({
    test_results$test5 <- run_test_5()
  }, error = function(e) {
    message("TEST 5 FAILED: ", e$message)
    test_results$test5 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  tryCatch({
    test_results$test6 <- run_test_6()
  }, error = function(e) {
    message("TEST 6 FAILED: ", e$message)
    test_results$test6 <- list(summary = list(success = FALSE, error = e$message))
  })
  
  # Generate summary report
  message("\n", paste(rep("=", 70), collapse=""))
  message("TEST SUITE SUMMARY")
  message(paste(rep("=", 70), collapse=""))
  
  success_count <- 0
  for (test_name in names(test_results)) {
    result <- test_results[[test_name]]
    status <- ifelse(result$summary$success, "PASS", "FAIL")
    message(sprintf("%-20s: %s", test_name, status))
    if (result$summary$success) success_count <- success_count + 1
  }
  
  message("\nTotal: ", success_count, "/", length(test_results), " tests passed")
  
  # Save summary
  summary_file <- file.path("/Users/bbuchsbaum/code/afnireg/afni_test", 
                            paste0("test_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
  saveRDS(test_results, summary_file)
  message("\nTest results saved to: ", summary_file)
  
  return(test_results)
}

# Function to run a single test by number
run_single_test <- function(test_num) {
  test_func <- switch(test_num,
    "1" = run_test_1,
    "2" = run_test_2,
    "3" = run_test_3,
    "4" = run_test_4,
    "5" = run_test_5,
    "6" = run_test_6,
    stop("Invalid test number. Use 1-6")
  )
  
  test_func()
}

# Print available tests
print_test_menu <- function() {
  cat("\nAvailable AFNI Tests:\n")
  cat("1. Simple single stimulus\n")
  cat("2. Two conditions (A/B)\n")
  cat("3. Multiple runs\n")
  cat("4. Different HRF models (BLOCK)\n")
  cat("5. Events with durations\n")
  cat("6. 2x2 Factorial design with contrasts\n")
  cat("\nUsage:\n")
  cat("  run_single_test(1)  # Run test 1\n")
  cat("  run_all_tests()     # Run all tests\n")
  cat("  view_test_results('path/to/test/dir')  # View results\n")
}

# Print menu on source
print_test_menu()