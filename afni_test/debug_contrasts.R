#!/usr/bin/env Rscript

# Debug script to understand how contrasts are passed to afni_hrf

library(fmrireg)
library(fmridesign)
library(fmrihrf)
library(afnireg)

# Set up test data
design <- data.frame(
  Fac1 = factor(c("aa", "bb", "aa", "bb")),
  Fac2 = factor(c("cc", "dd", "cc", "dd")),
  onset = c(10, 30, 50, 70),
  run = 1
)

con_fac1 <- pair_contrast(
  ~ Fac1 == "aa",
  ~ Fac1 == "bb",
  name = "Main_Fac1_aa_vs_bb"
)

all_contrasts <- contrast_set(con_fac1)

# Test what afni_hrf receives when called directly vs in formula
cat("=== DIRECT CALL TEST ===\n")
hrf_direct <- afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts)
cat("Direct call - contrasts in hrfspec:\n")
print(hrf_direct$contrasts)
print(class(hrf_direct$contrasts))

# Now test what happens when we evaluate the call using eval_tidy
cat("\n=== EVAL_TIDY TEST ===\n")
library(rlang)

# Create the expression as it would appear in the formula
expr <- quote(afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts))

# Create evaluation environment like find_and_eval_hrf_calls does
f_env <- environment()  # Current environment where all_contrasts is defined
eval_env <- rlang::env_bury(f_env, !!!design)

cat("Environment contents:\n")
cat("all_contrasts exists in f_env:", exists("all_contrasts", envir = f_env), "\n")
cat("all_contrasts exists in eval_env:", exists("all_contrasts", envir = eval_env), "\n")

# Evaluate the expression
hrf_eval <- rlang::eval_tidy(expr, env = eval_env)
cat("Eval_tidy call - contrasts in hrfspec:\n")
print(hrf_eval$contrasts)
print(class(hrf_eval$contrasts))

# Let's trace what afni_hrf receives for the contrasts parameter
cat("\n=== TRACING AFNI_HRF ARGUMENTS ===\n")

# Modify afni_hrf temporarily to see what it receives
original_afni_hrf <- afni_hrf

debug_afni_hrf <- function(..., basis=c("spmg1", "block", "dmblock",
                                  "tent",   "csplin", "poly",  "sin", "sine",
                                  "gam", "gamma", "spmg2", "spmg3", "wav"),
                                  onsets=NULL, durations=NULL, prefix=NULL, subset=NULL,
                                  nbasis=1, contrasts=NULL, id=NULL,
                                  lag=0, precision = 0.3, summate = TRUE,
                                  start=NULL, stop=NULL) {
  
  cat("DEBUG: contrasts parameter received:\n")
  print(contrasts)
  cat("DEBUG: contrasts class:\n")
  print(class(contrasts))
  cat("DEBUG: contrasts is.null:", is.null(contrasts), "\n")
  
  # Call original function
  result <- original_afni_hrf(..., basis=basis, onsets=onsets, durations=durations, 
                             prefix=prefix, subset=subset, nbasis=nbasis, 
                             contrasts=contrasts, id=id, lag=lag, precision=precision, 
                             summate=summate, start=start, stop=stop)
  
  cat("DEBUG: result contrasts:\n")
  print(result$contrasts)
  
  return(result)
}

# Replace afni_hrf temporarily
assign("afni_hrf", debug_afni_hrf, envir = .GlobalEnv)

cat("Testing with debug version:\n")
hrf_debug <- rlang::eval_tidy(expr, env = eval_env)

# Restore original function
assign("afni_hrf", original_afni_hrf, envir = .GlobalEnv)