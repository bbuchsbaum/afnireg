#!/usr/bin/env Rscript

# Simple test
library(testthat)
library(fmridesign)
library(fmrihrf)
source('tests/testthat/helper_afni.R')
library(afnireg)

design <- make_simple_design(n_cond = 2, n_runs = 1)
sframe <- make_sampling_frame(n_runs = 1)

model <- build_test_afni_model(
  design,
  onset ~ afni_hrf(condition, basis = 'SPMG1'),
  sframe
)

# Test the working approach
alm <- afni_lm(model, model$dataset, nodata = c(100, 2), x1D_stop = TRUE)
cat('Command contains -x1D_stop:', grepl('-x1D_stop', alm$cmd$cmd), '\n')

if (has_3ddeconvolve()) {
  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, 'afni_test_simple')
  
  result <- tryCatch({
    run.afni_lm_spec(alm, outdir = test_dir, execute = TRUE)
    TRUE
  }, error = function(e) {
    cat('Error:', e$message, '\n')
    FALSE
  })
  
  if (result && dir.exists(test_dir)) {
    files <- list.files(test_dir)
    cat('Files created:', length(files), '\n')
    print(files)
    
    # Look for any .1D files
    oneDfiles <- files[grepl('\\.1D$', files)]
    cat('1D files found:', length(oneDfiles), '\n')
    if (length(oneDfiles) > 0) {
      print(oneDfiles)
    }
  }
} else {
  cat('3dDeconvolve not available\n')
}