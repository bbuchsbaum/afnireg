#!/usr/bin/env Rscript

# Test script to investigate contrast evaluation in afni_hrf() 

library(fmrireg)
library(fmridesign)
library(fmrihrf)
library(afnireg)

# Create test data
design <- data.frame(
  Fac1 = factor(c("aa", "bb", "aa", "bb")),
  Fac2 = factor(c("cc", "dd", "cc", "dd")),
  onset = c(10, 30, 50, 70),
  run = 1
)

# Create contrasts
con_fac1 <- pair_contrast(
  ~ Fac1 == "aa",
  ~ Fac1 == "bb",
  name = "Main_Fac1_aa_vs_bb"
)

all_contrasts <- contrast_set(con_fac1)

print("Created all_contrasts:")
print(all_contrasts)
print("Class of all_contrasts:")
print(class(all_contrasts))

# Test 1: Direct passing of contrast object
cat("\n=== TEST 1: Direct contrast object ===\n")
try({
  hrf1 <- afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts)
  cat("SUCCESS: Direct contrast object worked\n")
  cat("Contrasts stored in hrf1:\n")
  print(hrf1$contrasts)
}, silent = FALSE)

# Test 2: Using variable name as symbol (like in formula)
cat("\n=== TEST 2: Variable name as symbol ===\n")
try({
  # This simulates what happens in formula evaluation
  # where 'all_contrasts' is a symbol that needs to be evaluated
  hrf2 <- afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = quote(all_contrasts))
  cat("SUCCESS: Variable name as symbol worked\n")
  cat("Contrasts stored in hrf2:\n")
  print(hrf2$contrasts)
}, silent = FALSE)

# Test 3: Create a formula with contrasts
cat("\n=== TEST 3: Formula with contrasts ===\n")
try({
  # This mimics what would happen in the formula
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  
  # Create the event model with the contrast variable in scope
  emodel <- event_model(onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts), 
                       data = design, 
                       block = ~ run, 
                       sampling_frame = sframe)
  
  cat("SUCCESS: Formula with contrasts worked\n")
  cat("Terms in event model:\n")
  print(names(terms(emodel)))
  
  # Check if contrasts made it through
  cat("Contrasts in first term:\n")
  first_term <- terms(emodel)[[1]]
  hrfspec <- attr(first_term, "hrfspec")
  if (!is.null(hrfspec) && !is.null(hrfspec$contrasts)) {
    print(hrfspec$contrasts)
  } else {
    cat("No contrasts found in hrfspec\n")
  }
  
}, silent = FALSE)

cat("\n=== TEST 4: Formula environment test ===\n")
# Test what happens when all_contrasts is in a different environment
test_formula_env <- function() {
  # all_contrasts is not defined in this function scope
  try({
    emodel <- event_model(onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts), 
                         data = design, 
                         block = ~ run, 
                         sampling_frame = sframe)
    cat("SUCCESS: Formula worked even without all_contrasts in local scope\n")
  }, error = function(e) {
    cat("ERROR: Formula failed when all_contrasts not in local scope:\n")
    cat(e$message, "\n")
  })
}

test_formula_env()