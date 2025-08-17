#' afnireg: AFNI Integration for fMRI Analysis in R
#'
#' The afnireg package provides an R interface to AFNI's 3dDeconvolve
#' program for fMRI analysis.
#'
#' @docType package
#' @name afnireg
#' @import fmridesign
#' @import fmrihrf
#' @import fmridataset
#' @import methods
#' @importFrom fmrireg fmri_model
#' @importFrom assertthat assert_that
#' @importFrom rlang enquo quo_get_expr
#' @importFrom purrr imap walk map map_chr
#' @importFrom utils read.table write.table
#' @importFrom stats terms formula model.matrix
NULL

# Re-export key functions from dependencies that are used in AFNI code
#' @export
fmridesign::event_model

#' @export
fmridesign::baseline_model

#' @export
fmridesign::event_term

#' @export
fmridesign::construct

#' @export
fmridesign::event_table

#' @export
fmridesign::split_onsets

#' @export
fmrihrf::blocklens

#' @export
fmridesign::blockids

#' @export
fmrihrf::sampling_frame

#' @export
fmridataset::fmri_dataset

#' @export
fmrireg::fmri_model
