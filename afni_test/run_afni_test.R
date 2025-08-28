#' Run AFNI Test with Structured Output
#' 
#' Executes an AFNI test case and captures all relevant information
#' in a structured format for inspection
#' 
#' @param test_name Character string identifying the test
#' @param design Data frame with event table
#' @param model_formula Formula for event model
#' @param sampling_frame Sampling frame object
#' @param baseline_model Baseline model object
#' @param test_description Character string describing the test
#' @param nodata Numeric vector c(NT, TR) for -nodata mode
#' @param contrasts Optional list of contrasts
#' @param options Additional options for afni_lm
#' @return List with test results and report
run_afni_test <- function(test_name, 
                         design, 
                         model_formula,
                         sampling_frame,
                         baseline_model,
                         test_description = "",
                         nodata = c(100, 2),
                         contrasts = NULL,
                         options = list()) {
  
  library(fmrireg)
  library(fmridesign)
  library(fmrihrf)
  library(afnireg)
  library(jsonlite)
  
  # Create test directory with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  test_dir <- file.path("/Users/bbuchsbaum/code/afnireg/afni_test", 
                        paste0(test_name, "_", timestamp))
  dir.create(test_dir, recursive = TRUE)
  
  # Store initial test specification
  test_spec <- list(
    test_name = test_name,
    description = test_description,
    timestamp = timestamp,
    date = Sys.time(),
    nodata = nodata,
    model_formula = deparse(model_formula),
    design_summary = list(
      n_events = nrow(design),
      variables = names(design),
      runs = if (is.null(design$run)) 1 else unique(design$run)
    )
  )
  
  # Initialize results
  results <- list(
    success = FALSE,
    errors = character(),
    warnings = character()
  )
  
  # Capture warnings
  warn_handler <- function(w) {
    results$warnings <<- c(results$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  
  # Build the model with warning capture
  withCallingHandlers({
    
    # Create event model
    emodel <- event_model(model_formula, 
                         data = design, 
                         block = ~ run, 
                         sampling_frame = sampling_frame)
    
    # Create dataset (dummy mode)
    dset <- fmri_dataset(
      scans = "dummy.nii", 
      mask = "dummy_mask.nii",
      TR = nodata[2],
      run_length = sampling_frame$blocklens[1],
      event_table = design,
      base_path = ".",
      dummy_mode = TRUE
    )
    
    # Create fmri_model
    fmodel <- fmri_model(emodel, baseline_model, dset)
    
    # Note: Contrasts should be provided in the model formula via hrf(contrasts=...)
    # not added separately here
    
    # Generate AFNI command
    message("Generating AFNI command...")
    alm <- afni_lm(fmodel, dset, nodata = nodata, options = options)
    
    # Store model components
    test_spec$model_components <- list(
      n_functional_terms = length(terms(emodel)),
      functional_term_names = names(terms(emodel)),
      baseline_terms = names(terms(baseline_model)),
      n_contrasts = length(contrasts)
    )
    
    # Store the command
    test_spec$afni_command <- alm$cmd$cmd
    test_spec$n_afni_stims <- length(alm$cmd$afni_stims)
    test_spec$n_baseline_mats <- length(alm$cmd$afni_baseline_mats)
    
    # Change to test directory for execution
    setwd(test_dir)
    
    # Write necessary files
    message("Writing stimulus files...")
    
    # Use the library's write functions which handle the format correctly
    afnireg:::write_stim_files(alm$cmd$afni_stims)
    
    # Write baseline matrices using the library function
    if (!is.null(alm$cmd$afni_baseline_mats)) {
      afnireg:::write_baseline_mats(alm$cmd$afni_baseline_mats)
    }
    
    # Write GLT files if any
    if (length(alm$cmd$gltfiles) > 0) {
      afnireg:::write_glts(alm$cmd$gltstr, alm$cmd$gltfiles)
    }
    
    # List files before execution
    files_before <- list.files()
    
    # Execute 3dDeconvolve
    message("Running 3dDeconvolve...")
    output <- system2("3dDeconvolve", 
                     args = strsplit(alm$cmd$cmd, " ")[[1]][-1],
                     stdout = TRUE, 
                     stderr = TRUE)
    
    # Store output
    test_spec$execution_output <- output
    
    # Parse output for key information
    test_spec$matrix_conditions <- list()
    test_spec$stimuli <- list()
    
    for (line in output) {
      # Check for matrix conditions
      if (grepl("matrix condition", line)) {
        condition <- gsub(".*\\[X\\][^:]*:\\s*([0-9.]+).*", "\\1", line)
        quality <- ifelse(grepl("VERY GOOD", line), "VERY GOOD",
                         ifelse(grepl("GOOD", line), "GOOD", "POOR"))
        test_spec$matrix_conditions[[length(test_spec$matrix_conditions) + 1]] <- 
          list(condition = condition, quality = quality)
      }
      
      # Check for stimulus info
      if (grepl("^Stimulus:", line)) {
        stim_name <- gsub("Stimulus:\\s*(.*)\\s*$", "\\1", line)
        test_spec$stimuli[[stim_name]] <- list(name = stim_name)
      }
      
      # Check for norm std dev
      if (grepl("norm\\. std\\. dev\\.", line)) {
        std_dev <- as.numeric(gsub(".*=\\s*([0-9.]+).*", "\\1", line))
        if (length(test_spec$stimuli) > 0) {
          test_spec$stimuli[[length(test_spec$stimuli)]]$std_dev <- std_dev
        }
      }
      
      # Check for errors (but not "average error" which is a metric)
      if (grepl("ERROR|FATAL", line, ignore.case = TRUE) && 
          !grepl("average error", line, ignore.case = TRUE)) {
        results$errors <- c(results$errors, line)
      }
    }
    
    # List files after execution
    files_after <- list.files()
    test_spec$files_created <- setdiff(files_after, files_before)
    
    # Check for key output files
    test_spec$xmat_file <- any(grepl("xmat\\.1D", files_after))
    test_spec$stat_files <- files_after[grep("statout|coefout", files_after)]
    
    # Determine success
    results$success <- length(results$errors) == 0 && 
                      test_spec$xmat_file &&
                      length(test_spec$matrix_conditions) > 0
    
    # Save model objects
    saveRDS(list(emodel = emodel, bmodel = baseline_model, 
                fmodel = fmodel, alm = alm), 
           file.path(test_dir, "model_objects.rds"))
    
  }, warning = warn_handler)
  
  # Save test specification
  test_spec$results <- results
  write_json(test_spec, file.path(test_dir, "test_spec.json"), 
            pretty = TRUE, auto_unbox = TRUE)
  
  # Write command to file
  writeLines(alm$cmd$cmd, file.path(test_dir, "command.txt"))
  
  # Write output to log
  writeLines(output, file.path(test_dir, "output.log"))
  
  # Generate summary
  summary <- list(
    test_name = test_name,
    timestamp = timestamp,
    success = results$success,
    n_warnings = length(results$warnings),
    n_errors = length(results$errors),
    xmat_created = test_spec$xmat_file,
    n_stimuli = length(test_spec$stimuli),
    test_dir = test_dir
  )
  
  message("\n=== Test Summary ===")
  message("Test: ", test_name)
  message("Success: ", results$success)
  message("Warnings: ", length(results$warnings))
  message("Errors: ", length(results$errors))
  message("Output directory: ", test_dir)
  
  return(list(
    summary = summary,
    test_spec = test_spec,
    test_dir = test_dir
  ))
}

# Convenience function to view test results
view_test_results <- function(test_dir) {
  spec <- fromJSON(file.path(test_dir, "test_spec.json"))
  
  cat("\n=== TEST RESULTS ===\n")
  cat("Test:", spec$test_name, "\n")
  cat("Description:", spec$description, "\n")
  cat("Date:", spec$date, "\n")
  cat("Success:", spec$results$success, "\n\n")
  
  cat("=== MODEL SPECIFICATION ===\n")
  cat("Formula:", spec$model_formula, "\n")
  cat("Functional terms:", paste(spec$model_components$functional_term_names, collapse=", "), "\n")
  cat("Baseline terms:", paste(spec$model_components$baseline_terms, collapse=", "), "\n\n")
  
  cat("=== AFNI COMMAND ===\n")
  cat(spec$afni_command, "\n\n")
  
  cat("=== MATRIX CONDITIONS ===\n")
  if (!is.null(spec$matrix_conditions) && length(spec$matrix_conditions) > 0) {
    for (cond in spec$matrix_conditions) {
      if (is.list(cond)) {
        cat(sprintf("Condition: %s - %s\n", cond$condition, cond$quality))
      }
    }
  }
  
  cat("\n=== STIMULI ===\n")
  for (stim in spec$stimuli) {
    if (!is.null(stim$std_dev)) {
      cat(sprintf("%s: std dev = %.4f\n", stim$name, stim$std_dev))
    } else {
      cat(sprintf("%s\n", stim$name))
    }
  }
  
  if (length(spec$results$warnings) > 0) {
    cat("\n=== WARNINGS ===\n")
    for (w in spec$results$warnings) {
      cat("- ", w, "\n")
    }
  }
  
  if (length(spec$results$errors) > 0) {
    cat("\n=== ERRORS ===\n")
    for (e in spec$results$errors) {
      cat("- ", e, "\n")
    }
  }
  
  cat("\n=== FILES CREATED ===\n")
  for (f in spec$files_created) {
    cat("- ", f, "\n")
  }
}