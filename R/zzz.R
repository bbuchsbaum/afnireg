#' Package Loading Hook
#'
#' Register AFNI HRF specifications with fmridesign when the package is loaded.
#'

.onLoad <- function(libname, pkgname) {
  # Check if fmridesign is available
  if (requireNamespace("fmridesign", quietly = TRUE)) {
    # Register AFNI HRF specifications with fmridesign
    if (exists("register_hrfspec_extension", where = asNamespace("fmridesign"))) {
      # Register afni_hrf specification
      fmridesign::register_hrfspec_extension(
        spec_class = "afni_hrfspec",
        package = "afnireg",
        convolved_class = "afni_hrf_convolved_term",
        requires_external_processing = TRUE
      )
      
      # Also register the convolved term class
      fmridesign::register_hrfspec_extension(
        spec_class = "afni_hrf_convolved_term",
        package = "afnireg",
        convolved_class = NULL,
        requires_external_processing = TRUE
      )
      
      # Register afni_trialwise specification
      fmridesign::register_hrfspec_extension(
        spec_class = "afni_trialwise_hrfspec",
        package = "afnireg",
        convolved_class = "afni_trialwise_convolved_term",
        requires_external_processing = TRUE
      )
      
      # Also register the trialwise convolved term class  
      fmridesign::register_hrfspec_extension(
        spec_class = "afni_trialwise_convolved_term",
        package = "afnireg",
        convolved_class = NULL,
        requires_external_processing = TRUE
      )
      
      packageStartupMessage("AFNI HRF specifications registered with fmridesign")
    }
  }
  
  invisible()
}

.onUnload <- function(libpath) {
  # Clean up could go here if needed
  invisible()
}