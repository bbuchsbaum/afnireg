#!/usr/bin/env Rscript

# Trace the options bug systematically

library(testthat)
library(fmridesign)
library(fmrihrf)
source('tests/testthat/helper_afni.R')
library(afnireg)

# Create a debug version of gen_afni_lm.fmri_model
debug_gen_afni_lm.fmri_model <- function(x, dataset, working_dir=".", 
                                    polort=-1, 
                                    jobs=1, 
                                    censor=NULL, 
                                    iresp=FALSE, 
                                    nodata=NULL,
                                    TR_times=x$sampling_frame$TR,
                                    x1D_stop=FALSE,
                                    nofullf_atall=TRUE,
                                    ...) {
  
  cat("=== debug_gen_afni_lm.fmri_model called ===\n")
  cat("x1D_stop parameter:", x1D_stop, "\n")
  
  # Check what's in ...
  dots <- list(...)
  cat("Arguments in ... :\n")
  print(str(dots))
  
  options <- list(...)
  cat("options after list(...):\n")
  print(str(options))
  
  # Set defaults for common options, but preserve any existing values in options
  if (is.null(options$polort)) options$polort <- polort
  if (is.null(options$jobs)) options$jobs <- jobs
  if (is.null(options$censor)) options$censor <- censor
  if (is.null(options$iresp)) options$iresp <- iresp
  if (is.null(options$nodata)) options$nodata <- nodata
  if (is.null(options$TR_times)) options$TR_times <- TR_times
  if (is.null(options$x1D_stop)) options$x1D_stop <- x1D_stop
  if (is.null(options$nofullf_atall)) options$nofullf_atall <- nofullf_atall
  
  cat("options after setting defaults:\n")
  print(str(options))
  cat("options$x1D_stop:", options$x1D_stop, "\n")
  
  # Skip the actual command building for now
  cmd <- list(cmd = "mock command", cmdlines = list(x1D_stop = options$x1D_stop))
  
  ret <- list(
    model = x,
    dataset = dataset,
    working_dir = working_dir,
    options = options,
    cmd = cmd
  )
  
  class(ret) <- "afni_lm_spec"
  ret
}

# Temporarily replace the function
assign("gen_afni_lm.fmri_model", debug_gen_afni_lm.fmri_model, envir = .GlobalEnv)

# Create test data
design <- make_simple_design(n_cond = 2, n_runs = 1)
sframe <- make_sampling_frame(n_runs = 1)

model <- build_test_afni_model(
  design,
  onset ~ afni_hrf(condition, basis = 'SPMG1'),
  sframe
)

cat("=== Testing with options list ===\n")
alm <- afni_lm(model, model$dataset, nodata = c(100, 2), 
               options = list(x1D_stop = TRUE))

cat("\n=== Final result ===\n")
print(str(alm$options))