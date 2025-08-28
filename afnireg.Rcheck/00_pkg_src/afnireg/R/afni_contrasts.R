#' AFNI Contrast Functions
#'
#' This file contains functions for handling contrasts in AFNI models,
#' including conversion to GLT (General Linear Test) format.
#'

# Note: The generic to_glt is defined in afni_generics.R


#' @export
to_glt.contrast <- function(x,...) {
  if (is.matrix(x$weights) && ncol(x$weights) > 1) {
    glts <- lapply(1:ncol(x$weights), function(i) {
      # Only include non-zero weights in GLT string
      non_zero <- which(x$weights[,i] != 0)
      if (length(non_zero) == 0) return("")
      
      # Format each term with appropriate sign
      terms <- sapply(non_zero, function(j) {
        w <- signif(x$weights[j,i], 4)
        if (j == non_zero[1]) {
          # First term: no plus sign if positive
          paste0(w, "*", x$condnames[j])
        } else if (w >= 0) {
          # Subsequent positive terms need plus sign
          paste0("+", w, "*", x$condnames[j])
        } else {
          # Negative terms already have minus sign
          paste0(w, "*", x$condnames[j])
        }
      })
      paste(terms, collapse = " ")
    })
    
    ret <- list(glt_str=glts,
                name=paste0("GLT_", x$name, "_", 1:ncol(x$weights)),
                con=x)
    
    class(ret) <- "glt_contrast_list"
    ret
  } else {
    # Only include non-zero weights in GLT string
    non_zero <- which(x$weights != 0)
    if (length(non_zero) == 0) {
      glt <- ""
    } else {
      # Format each term with appropriate sign
      terms <- sapply(non_zero, function(j) {
        w <- signif(x$weights[j], 4)
        if (j == non_zero[1]) {
          # First term: no plus sign if positive
          paste0(w, "*", x$condnames[j])
        } else if (w >= 0) {
          # Subsequent positive terms need plus sign
          paste0("+", w, "*", x$condnames[j])
        } else {
          # Negative terms already have minus sign
          paste0(w, "*", x$condnames[j])
        }
      })
      glt <- paste(terms, collapse = " ")
    }
    
    ret <- list(glt_str=glt,
                name=x$name,  # Don't add GLT_ prefix - tests expect original name
                con=x)
    
    class(ret) <- "glt"  # Use class "glt" not "glt_contrast"
    ret
  }
}

#' Extract contrasts from an afni_hrf_convolved_term
#' @keywords internal
#' @export
contrasts.afni_hrf_convolved_term <- function(x, ...) {
  # First check if contrasts are stored directly in the hrfspec
  if (!is.null(x$hrfspec) && !is.null(x$hrfspec$contrasts)) {
    return(x$hrfspec$contrasts)
  }
  
  # Otherwise check for contrasts attribute
  cset <- attr(x, "contrasts")
  if (!is.null(cset)) {
    return(cset)
  }
  
  # Return NULL if no contrasts found
  NULL
}

#' Extract contrast weights for an afni_hrf_convolved_term
#' @keywords internal
#' @export
#' @importFrom fmridesign contrast_weights
contrast_weights.afni_hrf_convolved_term <- function(x, ...) {
  # Get the contrast set
  cset <- contrasts(x)
  if (is.null(cset)) {
    return(NULL)
  }
  
  # Process each contrast by delegating to fmridesign's contrast_weights
  result <- list()
  
  for (i in seq_along(cset)) {
    con <- cset[[i]]
    
    # Use fmridesign's contrast_weights method for the specific contrast type
    # This handles all contrast types correctly (formula, unit, oneway, etc.)
    contrast_result <- tryCatch({
      fmridesign::contrast_weights(con, x)
    }, error = function(e) {
      warning(paste("Failed to compute contrast weights for", 
                   con$name %||% paste0("contrast_", i), 
                   ":", e$message))
      NULL
    })
    
    if (!is.null(contrast_result)) {
      # Store with the contrast name as key
      result[[contrast_result$name]] <- contrast_result
    }
  }
  
  result
}