#!/usr/bin/env Rscript

# Debug the assignment in afni_hrf

library(fmrireg)
library(fmridesign)
library(fmrihrf)
library(afnireg)

con_fac1 <- pair_contrast(
  ~ Fac1 == "aa",
  ~ Fac1 == "bb",
  name = "Main_Fac1_aa_vs_bb"
)

all_contrasts <- contrast_set(con_fac1)

# Recreate the exact logic from afni_hrf
cat("=== Recreating afni_hrf logic ===\n")

# The parameters that would be passed to afni_hrf
contrasts <- all_contrasts
basis <- "spmg1"
nbasis <- 1
onsets <- NULL
durations <- NULL  
prefix <- NULL
lag <- 0
precision <- 0.3
summate <- TRUE

# The logic from afni_hrf
basis <- tolower(basis)
basis <- match.arg(basis, c("spmg1", "block", "dmblock", "tent", "csplin", "poly", "sin", "sine", "gam", "gamma", "spmg2", "spmg3", "wav"))

# Simulate the vars processing (this would come from ...)
vars <- list(quote(Fac1), quote(Fac2))
vars_quos <- lapply(vars, function(expr) rlang::new_quosure(expr, env = rlang::caller_env()))
var_labels <- sapply(vars, rlang::as_label)
varnames <- sapply(vars, function(v) {
  if (rlang::is_symbol(v)) as.character(v) else make.names(rlang::as_label(v))
})
term <- vars_quos
label <- paste0("afni_hrf(", paste0(var_labels, collapse=","), ")")

# The HRF creation
hrf <- afnireg:::get_AFNI_HRF(basis, nbasis=nbasis)

# The varnames logic
if (!is.null(prefix)) {
  varnames <- paste0(prefix, "_", varnames)
}

termname <- paste0(varnames, collapse="::")

id <- termname  # since id is NULL

# The contrasts logic
cat("Input contrasts:\n")
print(contrasts)
cat("Class:", class(contrasts), "\n")

cset <- if (!is.null(contrasts) && inherits(contrasts, "contrast_spec")) {
  cat("Taking branch 1: single contrast_spec\n")
  list(contrasts)
} else if (!is.null(contrasts) && inherits(contrasts, "contrast_set")) {
  cat("Taking branch 2: contrast_set\n")
  contrasts
} else if (!is.null(contrasts) && is.list(contrasts)) {
  cat("Taking branch 3: generic list\n")
  if (all(sapply(contrasts, inherits, "contrast_spec"))) {
    contrasts
  } else {
    NULL
  }
} else { 
  cat("Taking branch 4: default NULL\n")
  NULL 
}

cat("cset result:\n")
print(cset)
cat("cset class:", class(cset), "\n")
cat("cset is.null:", is.null(cset), "\n")

# Create the result list
ret <- list(
  name=termname,
  id=id,
  varnames=varnames,
  vars=term,
  label=label,
  hrf=hrf,
  onsets=onsets,
  durations=durations,
  prefix=prefix,
  subset=NULL, # subset would be NULL in this case
  lag=lag,
  precision = precision,
  summate = summate,
  contrasts=cset)

cat("\nFinal ret$contrasts:\n")
print(ret$contrasts)
cat("ret$contrasts class:", class(ret$contrasts), "\n")
cat("ret$contrasts is.null:", is.null(ret$contrasts), "\n")

# Set class
class(ret) <- c("afni_hrfspec", "hrfspec", "list")

cat("\nFinal ret object contrasts:\n")
print(ret$contrasts)