#' Construct Methods for AFNI HRF Specifications
#'
#' This file contains the construct methods that build AFNI convolved terms
#' from HRF specifications and model data.
#'
#' @importFrom fmridesign event_term construct
#' @importFrom rlang quo_get_expr

#' Construct an afni_hrf_convolved_term from an afni_hrfspec
#' 
#' @param x An afni_hrfspec object
#' @param model_spec Model specification containing data and sampling frame
#' @param ... Additional arguments
#' @export
#' @importFrom fmridesign construct
construct.afni_hrfspec <- function(x, model_spec, ...) {
  # Extract variables from the model data
  vars <- x$vars
  data <- model_spec$data
  
  # Check if we're dealing with a subset
  if (!is.null(x$subset)) {
    subset_idx <- eval(x$subset, envir = data, enclos = model_spec$formula_env)
    data <- data[subset_idx, , drop = FALSE]
    model_spec$onsets <- model_spec$onsets[subset_idx]
    model_spec$blockids <- model_spec$blockids[subset_idx]
    if (!is.null(model_spec$durations)) {
      model_spec$durations <- model_spec$durations[subset_idx]
    }
  }
  
  # Evaluate the variables to get the actual data
  var_data <- lapply(vars, function(v) {
    eval(rlang::quo_get_expr(v), envir = data, enclos = model_spec$formula_env)
  })
  
  # Fix: Add variable names from the hrfspec
  names(var_data) <- x$varnames
  
  # Create the event term
  evterm <- event_term(
    evlist = var_data,
    onsets = model_spec$onsets,
    blockids = model_spec$blockids,
    durations = model_spec$durations
  )
  
  # Create the convolved term
  ret <- list(
    evterm = evterm,
    hrfspec = x,
    varname = x$name,
    sampling_frame = model_spec$sampling_frame
  )
  
  class(ret) <- c("afni_hrf_convolved_term", "event_term")
  
  # Store the hrfspec as an attribute for later access
  attr(ret, "hrfspec") <- x
  
  ret
}

#' Construct an afni_trialwise_convolved_term from an afni_trialwise_hrfspec
#' 
#' @param x An afni_trialwise_hrfspec object
#' @param model_spec Model specification containing data and sampling frame
#' @param ... Additional arguments
#' @export
construct.afni_trialwise_hrfspec <- function(x, model_spec, ...) {
  # Extract the modulation variable
  data <- model_spec$data
  
  # Check if we're dealing with a subset
  if (!is.null(x$subset)) {
    subset_idx <- eval(x$subset, envir = data, enclos = model_spec$formula_env)
    data <- data[subset_idx, , drop = FALSE]
    model_spec$onsets <- model_spec$onsets[subset_idx]
    model_spec$blockids <- model_spec$blockids[subset_idx]
    if (!is.null(model_spec$durations)) {
      model_spec$durations <- model_spec$durations[subset_idx]
    }
  }
  
  # For trialwise, we need a modulation value
  # This should be passed in the data with the same name as the label
  varname <- x$varname
  
  if (varname %in% names(data)) {
    value <- data[[varname]]
  } else {
    # Default to all 1s if no modulation variable found
    value <- rep(1, nrow(data))
  }
  
  # Create the event term with modulation
  evterm <- event_term(
    evlist = list(value = value),
    onsets = model_spec$onsets,
    blockids = model_spec$blockids,
    durations = model_spec$durations
  )
  
  # Create the trialwise convolved term
  ret <- list(
    evterm = evterm,
    hrfspec = x,
    varname = x$varname,
    sampling_frame = model_spec$sampling_frame
  )
  
  class(ret) <- c("afni_trialwise_convolved_term", "event_term")
  
  # Store the hrfspec as an attribute
  attr(ret, "hrfspec") <- x
  
  ret
}