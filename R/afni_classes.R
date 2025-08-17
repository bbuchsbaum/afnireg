#' AFNI-specific Convolved Term Classes
#'
#' This file contains the AFNI-specific convolved term classes that were
#' previously defined in fmridesign. These classes represent terms that
#' are processed by AFNI's 3dDeconvolve rather than R's convolution.
#'

#' Print method for afni_hrf_convolved_term objects.
#'
#' @param x An afni_hrf_convolved_term object.
#' @param ... Additional arguments.
#' @export
#' @rdname print
print.afni_hrf_convolved_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  Conditions: ", conditions(x), "\n")
  cat("  Term Types: ", paste(purrr::map_chr(x$evterm$events, ~ class(.)[[1]])), "\n")
}

#' Print method for afni_trialwise_convolved_term objects.
#'
#' @param x An afni_trialwise_convolved_term object.
#' @param ... Additional arguments.
#' @export
#' @rdname print
print.afni_trialwise_convolved_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  Conditions: ", conditions(x), "\n")
  cat("  Term Types: ", paste(purrr::map_chr(x$evterm$events, ~ class(.)[[1]])), "\n")
}

#' Check if an Object is an AFNI Convolved Term
#'
#' @param x An object to check
#' @return Logical indicating if the object is an AFNI convolved term
#' @export
is_afni_convolved_term <- function(x) {
  inherits(x, "afni_hrf_convolved_term") || 
  inherits(x, "afni_trialwise_convolved_term")
}

#' Extract Conditions from AFNI Convolved Terms
#'
#' @param x An AFNI convolved term object
#' @param ... Additional arguments
#' @return Character vector of condition names
#' @export
conditions.afni_hrf_convolved_term <- function(x, ...) {
  if (!is.null(x$evterm)) {
    conditions(x$evterm, ...)
  } else {
    character(0)
  }
}

#' @export
conditions.afni_trialwise_convolved_term <- function(x, ...) {
  if (!is.null(x$evterm)) {
    conditions(x$evterm, ...)
  } else {
    character(0)
  }
}

#' Extract Design Matrix from AFNI Convolved Terms
#'
#' AFNI convolved terms don't produce R design matrices as they are
#' processed externally by 3dDeconvolve.
#'
#' @param x An AFNI convolved term object
#' @param ... Additional arguments
#' @return NULL or empty matrix
#' @export
design_matrix.afni_hrf_convolved_term <- function(x, ...) {
  # AFNI terms are processed externally, return empty matrix
  warning("AFNI terms don't produce R design matrices. Use AFNI's 3dDeconvolve for processing.", 
          call. = FALSE)
  return(NULL)
}

#' @export
design_matrix.afni_trialwise_convolved_term <- function(x, ...) {
  # AFNI terms are processed externally, return empty matrix
  warning("AFNI trialwise terms don't produce R design matrices. Use AFNI's 3dDeconvolve for processing.", 
          call. = FALSE)
  return(NULL)
}

#' Check if a Term Requires External Processing
#'
#' @param x An object to check
#' @return Logical indicating if external processing is required
#' @export
requires_external_processing <- function(x) {
  UseMethod("requires_external_processing")
}

#' @export
requires_external_processing.default <- function(x) {
  FALSE
}

#' @export
requires_external_processing.afni_hrf_convolved_term <- function(x) {
  TRUE
}

#' @export
requires_external_processing.afni_trialwise_convolved_term <- function(x) {
  TRUE
}

#' Get Number of Basis Functions for AFNI Terms
#'
#' @param x An AFNI convolved term object
#' @return Number of basis functions
#' @export
nbasis.afni_hrf_convolved_term <- function(x) {
  # Get the actual nbasis from the HRF object
  if (!is.null(x$hrfspec) && !is.null(x$hrfspec$hrf)) {
    # The AFNI_HRF object has nbasis as an attribute
    nbasis_val <- attr(x$hrfspec$hrf, "nbasis")
    if (!is.null(nbasis_val)) {
      return(as.integer(nbasis_val))
    }
  }
  # Default to 1 if no HRF or nbasis found
  1
}

#' @export
nbasis.afni_trialwise_convolved_term <- function(x) {
  # Get the actual nbasis from the HRF object
  if (!is.null(x$hrfspec) && !is.null(x$hrfspec$hrf)) {
    # The AFNI_HRF object has nbasis as an attribute
    nbasis_val <- attr(x$hrfspec$hrf, "nbasis")
    if (!is.null(nbasis_val)) {
      return(as.integer(nbasis_val))
    }
  }
  # Default to 1 if no HRF or nbasis found
  1
}

#' Get Short Names for AFNI Terms
#'
#' @param x An AFNI convolved term object
#' @return Character vector of short names
#' @export
shortnames.afni_hrf_convolved_term <- function(x) {
  # For AFNI terms, shortnames match longnames
  longnames(x)
}

#' @export
shortnames.afni_trialwise_convolved_term <- function(x) {
  # For AFNI terms, shortnames match longnames
  longnames(x)
}

#' Check if AFNI Term is Continuous
#'
#' @param x An AFNI convolved term object
#' @return Logical indicating if the term is continuous
#' @export
#' @importFrom fmridesign is_continuous
is_continuous.afni_hrf_convolved_term <- function(x) {
  # Check if the underlying variable is continuous
  if (!is.null(x$evterm)) {
    # Check if the event term has levels (categorical)
    if (!is.null(x$evterm$levels)) {
      return(FALSE)
    }
    # Check if the event term has multiple conditions
    conds <- conditions(x)
    if (length(conds) > 1) {
      return(FALSE)
    }
  }
  # Default to FALSE for categorical
  FALSE
}

#' @export
is_continuous.afni_trialwise_convolved_term <- function(x) {
  # Trialwise terms are typically continuous (individual trial modulation)
  TRUE
}

#' Get Number of Basis Functions for Event Terms
#' 
#' @param x An event_term object  
#' @return Number of basis functions
#' @export
#' @importFrom fmrihrf nbasis
nbasis.event_term <- function(x) {
  # For R-based event_term, check the hrfspec attribute
  hrfspec <- attr(x, "hrfspec")
  if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) {
    # If HRF has nbasis attribute, use it
    hrf_nbasis <- attr(hrfspec$hrf, "nbasis")
    if (!is.null(hrf_nbasis)) {
      return(as.integer(hrf_nbasis))
    }
  }
  # Default to 1 if no nbasis info found
  1
}