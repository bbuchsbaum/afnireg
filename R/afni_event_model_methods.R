#' Contrast weights for event models with AFNI terms
#'
#' This method extends the fmridesign contrast_weights method to handle
#' event models that may contain AFNI terms. It filters out AFNI terms
#' and only processes R-based terms.
#'
#' @param x An event_model object
#' @param ... Additional arguments passed to methods
#' @return A list of contrast weights from R-based terms only, or an empty list
#' @export
#' @importFrom fmridesign contrast_weights
#' @importFrom stats terms
contrast_weights.event_model <- function(x, ...) {
  # Get all terms
  model_terms <- terms(x)
  
  # Filter to only R-based terms (non-AFNI)
  r_terms <- model_terms[!sapply(model_terms, requires_external_processing)]
  
  if (length(r_terms) == 0) {
    # No R-based terms, return empty list
    return(list())
  }
  
  # Extract contrast weights from R-based terms
  weights_list <- lapply(r_terms, function(term) {
    tryCatch({
      contrast_weights(term)
    }, error = function(e) {
      NULL
    })
  })
  
  # Remove NULL entries and flatten
  weights_list <- weights_list[!sapply(weights_list, is.null)]
  
  if (length(weights_list) == 0) {
    return(list())
  }
  
  # Flatten the list if needed
  do.call(c, weights_list)
}