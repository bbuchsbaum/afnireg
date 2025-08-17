#' S3 Methods for AFNI Terms
#'
#' This file contains S3 method implementations for AFNI-specific term classes,
#' providing standard interfaces for longnames, shortnames, cells, and other operations.
#'
#' @importFrom fmridesign longnames shortnames conditions cells is_continuous

#' @keywords internal
#' @export
longnames.afni_hrf_convolved_term <- function(x, ...) {
  # Get the event term conditions
  if (!is.null(x$evterm)) {
    conds <- conditions(x$evterm)
    # Return conditions as-is, preserving the format from event_term
    return(conds)
  }
  
  character(0)
}

#' Get short names for afni_hrf_convolved_term
#' 
#' @param x An afni_hrf_convolved_term object
#' @param ... Additional arguments
#' @importFrom fmridesign shortnames
#' @export
shortnames.afni_hrf_convolved_term <- function(x, ...) {
  # For factorial designs, shortnames should match longnames
  # since each combination is a separate regressor in AFNI
  longnames(x, ...)
}

#' Check if afni_hrf_convolved_term is continuous
#' 
#' @param x An afni_hrf_convolved_term object
#' @param ... Additional arguments
#' @export
is_continuous.afni_hrf_convolved_term <- function(x, ...) {
  # Delegate to the underlying event term
  if (!is.null(x$evterm)) {
    is_continuous(x$evterm, ...)
  } else {
    FALSE  # Default to FALSE (categorical)
  }
}

#' Get cells for afni_hrf_convolved_term
#' 
#' @param x An afni_hrf_convolved_term object
#' @param ... Additional arguments
#' @export
cells.afni_hrf_convolved_term <- function(x, ...) {
  # Delegate to the underlying event term
  if (!is.null(x$evterm)) {
    cells(x$evterm, ...)
  } else {
    NULL
  }
}