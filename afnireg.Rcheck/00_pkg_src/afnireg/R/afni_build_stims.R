#' AFNI Build Stims Functions
#'
#' This file contains functions for building AFNI stimulus objects from various term types.
#' These functions convert fmridesign terms into AFNI-compatible stimulus specifications.
#'
#' @importFrom fmridesign design_matrix split_onsets blockids cells conditions longnames
#' @importFrom assertthat assert_that

#' @keywords internal
#' @noRd
build_baseline_stims <- function(x) {
  if (!is.null(x$baseline_model)) {
    dmat <- design_matrix(x$baseline_model)
    labels <- paste0("baseline_", 1:ncol(dmat))
    baseline_files <- paste0(labels, ".1D")
    
    lapply(1:ncol(dmat), function(i) {
      structure(
        list(label=labels[i], file_name=baseline_files[i], values=dmat[,i]),
        class=c("afni_baseline_matrix", "afni_stim")
      )
    })
  } else {
    list()
  }
}

# Note: The generic build_afni_stims is defined in afni_generics.R

#' @keywords internal
#' @export
#' @importFrom fmridesign design_matrix longnames
build_afni_stims.convolved_term <- function(x, iresp=FALSE, tr_times=1) {
  # For R-based convolved terms, generate stim_file entries from design matrix columns
  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
  
  # Get the design matrix for this term
  desmat <- tryCatch({
    design_matrix(x)
  }, error = function(e) {
    warning("Could not extract design matrix for term: ", e$message)
    return(NULL)
  })
  
  if (is.null(desmat) || ncol(desmat) == 0) {
    return(NULL)
  }
  
  # Create an afni_stim_file for each column
  lapply(1:length(stimlabels), function(i) {
    structure(
      list(label=stimlabels[i], file_name=stimfiles[i], values=desmat[, i]),
      class=c("afni_stim_file", "afni_stim")
    )
  })
}

#' @keywords internal
#' @export
build_afni_stims.afni_hrf_convolved_term <- function(x, iresp=FALSE, tr_times=1) {

  stimlabels <- longnames(x)
  stimfiles <- paste(stimlabels, "_times.1D", sep = "")
  # Note: We don't need design_matrix for AFNI terms as they're handled externally
  
  hrf_name <- as.character(x$hrfspec$hrf)
  blids <- unique(blockids(x$evterm))
  # For factorial designs, split_onsets might return fewer elements than longnames
  # We need to handle this carefully
  split_ons <- split_onsets(x$evterm, x$sampling_frame, global=FALSE, blocksplit = TRUE)
  
  # If the number of split onsets doesn't match stimlabels, 
  # we need to map them correctly
  if (length(split_ons) != length(stimlabels)) {
    # Get the actual condition names from the split
    actual_names <- names(split_ons)
    if (is.null(actual_names)) {
      # If no names, use the cells from evterm
      evterm_cells <- cells(x$evterm)
      if (!is.null(evterm_cells) && length(evterm_cells) == length(split_ons)) {
        actual_names <- evterm_cells
      }
    }
    
    # Create a properly sized list with all stimlabels
    new_split_ons <- vector("list", length(stimlabels))
    names(new_split_ons) <- stimlabels
    
    # Fill in the onsets we have
    for (i in seq_along(split_ons)) {
      if (!is.null(actual_names) && actual_names[i] %in% stimlabels) {
        new_split_ons[[actual_names[i]]] <- split_ons[[i]]
      } else if (i <= length(stimlabels)) {
        new_split_ons[[i]] <- split_ons[[i]]
      }
    }
    
    # For any missing conditions, use empty onset lists
    nruns <- length(unique(blockids(x$evterm)))
    for (i in seq_along(new_split_ons)) {
      if (is.null(new_split_ons[[i]])) {
        # Create empty onset structure with proper number of runs
        new_split_ons[[i]] <- rep(list(numeric(0)), nruns)
      } else if (length(new_split_ons[[i]]) < nruns) {
        # Pad with empty runs if condition doesn't appear in all runs
        while (length(new_split_ons[[i]]) < nruns) {
          new_split_ons[[i]] <- c(new_split_ons[[i]], list(numeric(0)))
        }
      }
    }
    
    split_ons <- new_split_ons
  } else {
    names(split_ons) <- stimlabels
    # Also need to pad here
    nruns <- length(unique(blockids(x$evterm)))
    for (i in seq_along(split_ons)) {
      if (length(split_ons[[i]]) < nruns) {
        while (length(split_ons[[i]]) < nruns) {
          split_ons[[i]] <- c(split_ons[[i]], list(numeric(0)))
        }
      }
    }
  }
  
  # Get durations if specified
  durations <- x$hrfspec$durations
  
  ret <- lapply(1:length(stimlabels), function(i) {
    stim_data <- list(label=stimlabels[i], file_name=stimfiles[i], hrf=hrf_name, 
                      onsets=split_ons[[stimlabels[[i]]]], onset_list=split_ons[[stimlabels[[i]]]], 
                      blockids=blids, iresp=iresp, sresp=FALSE, tr_times=tr_times)
    
    # Add durations if they exist
    if (!is.null(durations)) {
      # If durations is a single value, replicate for each run
      if (length(durations) == 1) {
        stim_data$durations <- rep(list(durations), length(split_ons[[stimlabels[[i]]]]))
      } else {
        # Otherwise use as is
        stim_data$durations <- durations
      }
    }
    
    structure(stim_data, class=c("afni_stim_times","afni_stim"))
  })
  
}

#' @keywords internal
#' @export
build_afni_stims.afni_trialwise_convolved_term <- function(x, iresp=FALSE, tr_times=1) {
  eterm <- x$evterm
  sf <- x$sampling_frame
  
  stimlabel <- x$varname
  stimfile <- paste(stimlabel, "_times.1D", sep = "")
  
  hrf_name <- as.character(x$hrfspec$hrf)
  
  # Get the modulation values
  if (is.factor(eterm$value)) {
    stop("trialwise afni_hrf does not support factors, convert to numeric or use standard 'afni_hrf'")
  }
  
  modvals <- as.numeric(eterm$value)
  
  ## scale to max abs == 1
  maxval <- max(abs(modvals))
  if (maxval > 0) {
    modvals <- modvals/maxval
  }
  
  # Split onsets and modulation values by run
  onsets_by_run <- split_onsets(eterm, sf, global=FALSE, blocksplit=TRUE)[[1]]
  mods_by_run <- split(modvals, blockids(eterm))
  
  # Format as "onset*modvalue" for each event
  formatted_onsets <- mapply(function(ons, mods) {
    if (length(ons) > 0) {
      paste(paste0(ons, "*", mods), collapse=" ")
    } else {
      "*"  # Empty run
    }
  }, onsets_by_run, mods_by_run, SIMPLIFY = FALSE)
  
  list(structure(
    list(label=stimlabel, file_name=stimfile, hrf=hrf_name, 
         onsets=formatted_onsets, blockids=unique(blockids(eterm))),
    class=c("afni_stim_im_times","afni_stim")
  ))
}

#' @keywords internal
#' @export
#' @importFrom fmridesign split_onsets
build_afni_stims.event_term <- function(x, iresp=FALSE, tr_times=1,...) {
  # Check if this is an R-based convolved term (has hrfspec but not AFNI)
  hrfspec <- attr(x, "hrfspec")
  if (!is.null(hrfspec) && !inherits(hrfspec, "afni_hrfspec") && !inherits(hrfspec, "afni_trialwise_hrfspec")) {
    # This is an R-based term, treat it like convolved_term
    # Generate stim_file entries from design matrix columns
    stimlabels <- longnames(x)
    stimfiles <- paste(stimlabels, "_reg.1D", sep = "")
    
    # Get the design matrix for this term
    desmat <- tryCatch({
      design_matrix(x)
    }, error = function(e) {
      warning("Could not extract design matrix for event_term: ", e$message)
      return(NULL)
    })
    
    if (is.null(desmat) || ncol(desmat) == 0) {
      return(NULL)
    }
    
    # Create an afni_stim_file for each column
    # Use the internal afni_stim_file function from afni_stim_files.R
    stims <- lapply(1:length(stimlabels), function(i) {
      structure(
        list(label=stimlabels[i], file_name=stimfiles[i], values=desmat[, i]),
        class=c("afni_stim_file", "afni_stim")
      )
    })
    return(stims)
  } else if (length(x$events) > 0) {
    # Otherwise handle compound event terms that might contain AFNI terms
    # Process each event in the compound term
    stims <- lapply(x$events, function(evt) {
      if (inherits(evt, "afni_hrf_convolved_term") || 
          inherits(evt, "afni_trialwise_convolved_term")) {
        build_afni_stims(evt, iresp=iresp, tr_times=tr_times)
      } else {
        NULL
      }
    })
    
    # Filter out NULLs and unlist one level
    stims <- Filter(Negate(is.null), stims)
    if (length(stims) > 0) {
      unlist(stims, recursive = FALSE)
    } else {
      NULL
    }
  } else {
    NULL
  }
}