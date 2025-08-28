#!/usr/bin/env Rscript

# Debug the inheritance chain for contrast_set

library(fmridesign)

con_fac1 <- pair_contrast(
  ~ Fac1 == "aa",
  ~ Fac1 == "bb",
  name = "Main_Fac1_aa_vs_bb"
)

all_contrasts <- contrast_set(con_fac1)

cat("Class hierarchy for all_contrasts:\n")
print(class(all_contrasts))

cat("\nChecking inheritance:\n")
cat("inherits(all_contrasts, 'contrast_spec'):", inherits(all_contrasts, "contrast_spec"), "\n")
cat("inherits(all_contrasts, 'contrast_set'):", inherits(all_contrasts, "contrast_set"), "\n")
cat("is.list(all_contrasts):", is.list(all_contrasts), "\n")

cat("\nContents of all_contrasts:\n")
print(str(all_contrasts))

cat("\nIndividual elements:\n")
for (i in seq_along(all_contrasts)) {
  cat("Element", i, "class:", class(all_contrasts[[i]]), "\n")
  cat("Element", i, "inherits contrast_spec:", inherits(all_contrasts[[i]], "contrast_spec"), "\n")
}

cat("\nTesting the afni_hrf logic manually:\n")
contrasts <- all_contrasts

cset <- if (!is.null(contrasts) && inherits(contrasts, "contrast_spec")) {
  cat("Branch 1: Single contrast_spec\n")
  list(contrasts)
} else if (!is.null(contrasts) && inherits(contrasts, "contrast_set")) {
  cat("Branch 2: contrast_set\n")
  contrasts
} else if (!is.null(contrasts) && is.list(contrasts)) {
  cat("Branch 3: list - checking if all elements are contrast_spec\n")
  all_are_contrast_spec <- all(sapply(contrasts, inherits, "contrast_spec"))
  cat("All elements are contrast_spec:", all_are_contrast_spec, "\n")
  if (all_are_contrast_spec) {
    contrasts
  } else {
    NULL
  }
} else { 
  cat("Branch 4: Default to NULL\n")
  NULL 
}

cat("\nResult of cset assignment:\n")
print(cset)
print(class(cset))