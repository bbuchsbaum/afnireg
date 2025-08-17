#' AFNI Stimulus File Functions
#'
#' This file contains functions for creating, writing, and managing AFNI stimulus files
#' including stim_file, stim_times, stim_times_IM, and baseline matrices.
#'

#' @keywords internal
#' @noRd
afni_stim_file <- function(label, file_name, values) {
  structure(
    list(label=label, file_name=file_name, values=values),
    class=c("afni_stim_file", "afni_stim")
  )
}

#' @keywords internal
#' @noRd
afni_stim_times <- function(label, file_name, hrf, onsets, blockids, iresp=FALSE, sresp=FALSE, tr_times=1) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, onset_list=onsets, blockids=blockids, iresp=iresp, sresp=sresp, tr_times=tr_times),
    class=c("afni_stim_times","afni_stim")
  )
}


#' @keywords internal
#' @noRd
afni_stim_im_times <- function(label, file_name, hrf, onsets, blockids) {
  structure(
    list(label=label, file_name=file_name, hrf=hrf, onsets=onsets, blockids=blockids),
    class=c("afni_stim_im_times","afni_stim")
  )
}


#' @keywords internal
#' @export
afni_command_switch <- function(x, ...) UseMethod("afni_command_switch")

#' @keywords internal
#' @export
write_afni_stim <- function(x, ...) UseMethod("write_afni_stim")


#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_file <- function(x, k, type) {
  if (type == "file") {
    paste0("-stim_file ", k, " ", x$file_name)
  } else if (type == "label") {
    paste0("-stim_label ", k, " ", x$label)
  } else if (type == "ortvec") {
    paste0("-ortvec ", x$file_name, " ", x$label)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_im_times <- function(x, k, type) {
  if (type == "times_IM") {
    paste0("-stim_times_IM ", k, " ", x$file_name, " '", x$hrf, "'")
  } else if (type == "label") {
    paste0("-stim_label ", k, " ", x$label)
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
afni_command_switch.afni_stim_times <- function(x, k, type) {
  if (type == "times") {
    paste0("-stim_times ", k, " ", x$file_name, " '", x$hrf, "'")
  } else if (type == "label") {
    paste0("-stim_label ", k, " ", x$label)
  } else if (type == "iresp" && x$iresp) {
    paste0("-iresp ", k, " iresp_", x$label, ".1D")
  } else if (type == "sresp" && x$sresp) {
    paste0("-sresp ", k, " sresp_", x$label, ".1D")
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
next_dir_name <- function(wd) {
  nd <- paste(wd, "+", sep="")
  if (!file.exists(nd)) {
    nd
  } else {
    Recall(nd)
  }
}

#' @keywords internal
#' @noRd
write_baseline_mat <- function(stim, dir) {
  write.table(stim$values, paste(dir, stim$file_name, sep="/"), row.names=FALSE, col.names=FALSE)
}

#' @keywords internal
#' @noRd
write_baseline_mats <- function(blist) {
  lapply(blist, function(x) write_baseline_mat(x, "."))
}

#' @keywords internal
#' @noRd
write_stim_files <- function(afni_stims) {
  sapply(afni_stims, function(stim) {
    write_afni_stim(stim, ".")
  })
}

#' @keywords internal
#' @export
write_afni_stim.afni_stim_file <- function(stim, dir) {
  if (any(is.na(stim$values))) {
    stim$values[is.na(stim$values)] <- 0
  }
  
  # Handle both vectors and matrices
  if (is.matrix(stim$values) || is.data.frame(stim$values)) {
    write.table(stim$values, paste(dir, stim$file_name, sep="/"), row.names=FALSE, col.names=FALSE)
  } else {
    # For vectors, write as a single column
    write.table(matrix(stim$values, ncol=1), paste(dir, stim$file_name, sep="/"), row.names=FALSE, col.names=FALSE)
  }
}

#' @keywords internal
#' @export
write_afni_stim.afni_stim_times <- function(stim, dir) {
  # Create lines for each run
  lines <- sapply(stim$onsets, function(ons) {
    if (length(ons) == 0) {
      "*"  # Empty run marked with asterisk
    } else {
      paste(ons, collapse=" ")
    }
  }, USE.NAMES = FALSE)
  
  # Ensure lines is character vector
  lines <- as.character(lines)
  
  # Write to file
  writeLines(lines, paste(dir, stim$file_name, sep="/"))
}

#' @keywords internal
#' @export
write_afni_stim.afni_stim_im_times <- function(stim, dir) {
  # Similar to stim_times but with modulation values
  lines <- sapply(stim$onsets, function(ons) {
    if (length(ons) == 0) {
      "*"
    } else {
      paste(ons, collapse=" ")
    }
  })
  
  writeLines(lines, paste(dir, stim$file_name, sep="/"))
}

#' @keywords internal
#' @noRd
write_censor_file <- function(dir, censor) {
  censor <- ifelse(censor, 1, 0)
  write.table(censor, paste(dir, "censor.1D", sep="/"), row.names=FALSE, col.names=FALSE)
}

#' @keywords internal
#' @noRd
#' @importFrom utils write.table
write_glts <- function(glts, gltfiles) {
  lapply(seq_along(glts), function(i) {
    write(glts[i], file=gltfiles[i])
  })
}

#' @keywords internal
#' @noRd
afni_baseline_matrix <- function(label, file_name, mat) {
  structure(
    list(label=label, file_name=file_name, values=mat),
    class=c("afni_stim_file", "afni_stim")
  )
}