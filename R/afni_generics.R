# AFNI-specific generic functions and S3 methods
# Extracted from fmrireg/R/all_generic.R

#' Generate AFNI Linear Model
#' 
#' Generic function to generate AFNI linear model specifications
#' 
#' @param x The input object
#' @param ... Additional arguments
#' @export
gen_afni_lm <- function(x, ...) {
  UseMethod("gen_afni_lm")
}

#' Build AFNI Stimulus Files
#' 
#' Generic function to build AFNI stimulus files
#' 
#' @param x The input object
#' @param ... Additional arguments
#' @export
build_afni_stims <- function(x, ...) {
  UseMethod("build_afni_stims")
}

#' Convert to AFNI GLT
#' 
#' Generic function to convert contrasts to AFNI GLT format
#' 
#' @param x The input object
#' @param ... Additional arguments
#' @export
to_glt <- function(x, ...) {
  UseMethod("to_glt")
}

# S3 method registrations for AFNI classes will be handled by roxygen2
