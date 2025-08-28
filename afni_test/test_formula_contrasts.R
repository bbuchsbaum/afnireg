#!/usr/bin/env Rscript

# Test contrasts in formula context

library(devtools)
load_all()
library(fmrireg)
library(fmridesign)
library(fmrihrf)

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

cat("=== Testing contrasts in event_model formula ===\n")

sframe <- sampling_frame(blocklens = 100, TR = 2)

# Test with contrasts in formula
emodel <- event_model(onset ~ afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts), 
                     data = design, 
                     block = ~ run, 
                     sampling_frame = sframe)

cat("SUCCESS: Event model created with contrasts\n")
cat("Terms in event model:\n")
print(names(terms(emodel)))

# Check if contrasts made it through
cat("Checking contrasts in first term:\n")
first_term <- terms(emodel)[[1]]
hrfspec <- attr(first_term, "hrfspec")
if (!is.null(hrfspec) && !is.null(hrfspec$contrasts)) {
  cat("Contrasts found in hrfspec:\n")
  print(hrfspec$contrasts)
} else {
  cat("No contrasts found in hrfspec\n")
}

# Test GLT generation
cat("\n=== Testing GLT generation ===\n")
bmodel <- baseline_model(basis = "constant", sframe = sframe)
fmodel <- fmri_model(emodel, bmodel)

# Try to get contrast weights
tryCatch({
  cwts <- contrast_weights(fmodel)
  cat("Contrast weights extracted successfully:\n")
  print(names(cwts))
}, error = function(e) {
  cat("Error extracting contrast weights:", e$message, "\n")
})