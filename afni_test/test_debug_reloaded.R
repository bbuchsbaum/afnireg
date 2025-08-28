#!/usr/bin/env Rscript

# Test with reloaded debug version of afni_hrf

library(devtools)
load_all()
library(fmrireg)
library(fmridesign)
library(fmrihrf)

con_fac1 <- pair_contrast(
  ~ Fac1 == "aa",
  ~ Fac1 == "bb",
  name = "Main_Fac1_aa_vs_bb"
)

all_contrasts <- contrast_set(con_fac1)

cat("=== Testing debug afni_hrf with reloaded package ===\n")
hrf_result <- afni_hrf(Fac1, Fac2, basis = "SPMG1", contrasts = all_contrasts)

cat("\nResult contrasts:\n")
print(hrf_result$contrasts)