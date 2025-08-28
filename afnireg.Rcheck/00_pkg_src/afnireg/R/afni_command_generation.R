#' AFNI Command Generation Functions
#'
#' This file contains functions for building 3dDeconvolve command strings
#' and managing all the components needed for AFNI command execution.
#'
#' @importFrom fmridesign contrast_weights terms blockids
#' @importFrom assertthat assert_that

#' Generate 3dDeconvolve command string
#' 
#' @keywords internal
#' @noRd
.make_decon_command_str <- function(cmdlines) {
  
  cmd <- character()
  
  if (!is.null(cmdlines$nodata)) {
    cmd <- paste(cmd, "-nodata", cmdlines$nodata)
  }
  
  if (!is.null(cmdlines$concat)) {
    cmd <- paste(cmd, "-concat", paste0("'", cmdlines$concat, "'"))
  }
  
  if (!is.null(cmdlines$input)) {
    cmd <- paste(cmd, "-input", cmdlines$input)
  }
  
  if (!is.null(cmdlines$mask)) {
    cmd <- paste(cmd, "-mask", cmdlines$mask)
  }
  
  if (!is.null(cmdlines$censor)) {
    cmd <- paste(cmd, "-censor", cmdlines$censor)
  }
  
  if (!is.null(cmdlines$polort) && cmdlines$polort >= 0) {
    cmd <- paste(cmd, "-polort", cmdlines$polort)
  } else {
    cmd <- paste(cmd, "-polort", -1)
  }
  
  if (!is.null(cmdlines$global_times) && cmdlines$global_times == TRUE) {
    cmd <- paste(cmd, "-global_times")
  } else {
    cmd <- paste(cmd, "-local_times")
  }
  
  if (!is.null(cmdlines$TR_times)) {
    cmd <- paste(cmd, "-TR_times", cmdlines$TR_times)
  }
  
  cmd <- paste(cmd, "-num_stimts", cmdlines$num_stimts)
  
  if (!is.null(cmdlines$stim_file)) {
    for (el in cmdlines$stim_file) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$stim_times)) {
    for (el in cmdlines$stim_times) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$stim_times_IM)) {
    for (el in cmdlines$stim_times_IM) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$stim_label)) {
    for (el in cmdlines$stim_label) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$ortvec)) {
    for (el in cmdlines$ortvec) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$iresp)) {
    for (el in cmdlines$iresp) {
      cmd <- paste(cmd, el)
    }
  }
  
  if (!is.null(cmdlines$num_glt) && cmdlines$num_glt > 0) {
    cmd <- paste(cmd, "-num_glt", cmdlines$num_glt)
    
    if (!is.null(cmdlines$gltsym)) {
      for (i in 1:length(cmdlines$gltsym)) {
        cmd <- paste(cmd, "-gltsym 'SYM:", cmdlines$gltsym[[i]], "'", sep="")
      }
    }
    
    if (!is.null(cmdlines$glt_label)) {
      for (el in cmdlines$glt_label) {
        cmd <- paste(cmd, "-glt_label", el)
      }
    }
  }
  
  if (!is.null(cmdlines$nofullf_atall) && cmdlines$nofullf_atall == TRUE) {
    cmd <- paste(cmd, "-nofullf_atall")
  }
  
  if (!is.null(cmdlines$rout) && cmdlines$rout == TRUE) {
    cmd <- paste(cmd, "-rout")
  }
  
  if (!is.null(cmdlines$tout) && cmdlines$tout == TRUE) {
    cmd <- paste(cmd, "-tout")
  }
  
  if (!is.null(cmdlines$fout) && cmdlines$fout == TRUE) {
    cmd <- paste(cmd, "-fout")
  }
  
  if (!is.null(cmdlines$bout) && cmdlines$bout == TRUE) {
    cmd <- paste(cmd, "-bout")
  }
  
  if (!is.null(cmdlines$noFDR) && cmdlines$noFDR == TRUE) {
    cmd <- paste(cmd, "-noFDR")
  }
  
  if (!is.null(cmdlines$nocond) && cmdlines$nocond == TRUE) {
    cmd <- paste(cmd, "-nocond")
  }
  
  if (!is.null(cmdlines$x1D_stop) && cmdlines$x1D_stop == TRUE) {
    # Determine the x1D filename based on bucket name
    bucket_name <- if (!is.null(cmdlines$bucket)) cmdlines$bucket else "stats"
    x1D_filename <- paste0(bucket_name, ".xmat.1D")
    cmd <- paste(cmd, "-x1D", x1D_filename, "-x1D_stop")
  }
  
  if (!is.null(cmdlines$cbucket)) {
    cmd <- paste(cmd, "-cbucket", cmdlines$cbucket)
  }
  
  if (!is.null(cmdlines$bucket)) {
    cmd <- paste(cmd, "-bucket", cmdlines$bucket)
  } else {
    cmd <- paste(cmd, "-bucket", "stats")
  }
  
  if (!is.null(cmdlines$errts)) {
    cmd <- paste(cmd, "-errts", cmdlines$errts)
  }
  
  if (!is.null(cmdlines$jobs) && cmdlines$jobs > 1) {
    cmd <- paste(cmd, "-jobs", cmdlines$jobs)
  }
  
  if (!is.null(cmdlines$float) && cmdlines$float == TRUE) {
    cmd <- paste(cmd, "-float")
  }
  
  paste("3dDeconvolve", cmd, collapse=" ")
}


#' Build 3dDeconvolve command
#' 
#' This function builds a 3dDeconvolve command for AFNI based on the provided model, 
#' dataset, working directory, and options.
#'
#' @param model The fmri_model object
#' @param dataset The dataset object
#' @param working_dir The working directory
#' @param opts Options list
#'
#' @return A list containing:
#'         - cmd: The 3dDeconvolve command string
#'         - cmdlines: The command lines for the 3dDeconvolve command
#'         - afni_stims: A list of AFNI stimulus objects
#'         - afni_baseline_mats: A list of AFNI baseline matrices
#'         - gltfiles: A list of GLT (general linear test) filenames
#'         - gltnames: A list of GLT names
#'         - glts: A list of GLT objects
#'         - gltstr: A list of GLT strings
#'         - censor: The censoring vector
#'
#' @importFrom fmridesign contrast_weights
#' @keywords internal
#' @noRd
build_decon_command <- function(model, dataset, working_dir, opts) {
  func_terms <- terms(model$event_model)
  message("number of functional terms: ", length(func_terms))
  
  # First, generate AFNI stims and filter out NULLs
  afni_stims <- lapply(func_terms, function(term) { build_afni_stims(term, iresp=opts[["iresp"]], tr_times=opts[["TR_times"]]) })
  afni_stims <- Filter(Negate(is.null), afni_stims) # Filter out NULLs
  afni_stims <- unlist(afni_stims, recursive = FALSE) # Unlist one level
  
  # Now get stimlabels only for terms that generated AFNI stims
  # We need to identify which terms generated stims
  terms_with_stims <- func_terms[!sapply(lapply(func_terms, function(term) { build_afni_stims(term, iresp=opts[["iresp"]], tr_times=opts[["TR_times"]]) }), is.null)]
  stimlabels <- unlist(lapply(terms_with_stims, longnames))
  
  ## all stims must be unique
  assert_that(length(unique(stimlabels)) == length(stimlabels))

  # Note: We can't assert stimlabels == conditions because conditions includes ALL terms,
  # but stimlabels only includes terms that generate AFNI stims
  # assert_that(length(stimlabels) == length(conditions(model$event_model)))
  
  ## extract all contrast matrices
  ## First try the standard way
  cons <- tryCatch({
    contrast_weights(model)
  }, error = function(e) {
    message("Standard contrast extraction failed: ", e$message)
    NULL
  })
  
  message("Standard extraction returned ", 
          if(is.null(cons)) "NULL" else paste(length(cons), "contrasts"))
  
  ## If standard extraction failed or returned NULL, try extracting from AFNI terms directly
  if (is.null(cons) || length(cons) == 0) {
    message("Trying AFNI term extraction...")
    cons <- list()
    # Extract contrasts from AFNI terms that have them in hrfspec
    for (term in func_terms) {
      if (inherits(term, "afni_hrf_convolved_term")) {
        hrfspec <- attr(term, "hrfspec")
        if (!is.null(hrfspec) && !is.null(hrfspec$contrasts)) {
          # Get contrast weights for this term
          term_cwts <- contrast_weights(term)
          if (!is.null(term_cwts)) {
            # Extract just the contrast objects without preserving list names
            # to avoid prefixing issues
            for (cw in term_cwts) {
              cons <- c(cons, list(cw))
            }
            message("Extracted ", length(term_cwts), " contrasts from AFNI term")
          }
        }
      }
    }
  }
  
  ## ensure contrasts are valid objects
  cons <- Filter(function(x) inherits(x, "contrast"), cons)
  assert_that(is.list(cons), all(sapply(cons, inherits, "contrast")))
  
  ## convert to 'glt's
  glts <- lapply(cons, to_glt)
  
  gltfiles <- unlist(lapply(glts, function(x) paste0(x$name, ".txt")), use.names = FALSE)
  gltnames <- unlist(lapply(glts, function(x) x$name), use.names = FALSE)
  gltstr <- unlist(lapply(glts, function(x) x$glt_str), use.names = FALSE)
  
  assert_that(sum(duplicated(gltnames))  == 0, msg="Cannot have two GLTs with the same name")
  
  afni_baseline_mats <- build_baseline_stims(model)
  
  purge_nulls <- function(A) {
    A[!sapply(A, is.null)]
  }
  
  opt_stim_labels <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "label")))
  opt_stim_files  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "file")))
  opt_stim_times  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "times")))
  opt_stim_times_IM  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "times_IM")))
  opt_stim_ortvecs <- purge_nulls(lapply(seq_along(afni_baseline_mats), function(i) afni_command_switch.afni_stim_file(afni_baseline_mats[[i]], i, "ortvec")))
  opt_stim_iresp  <-  purge_nulls(lapply(seq_along(afni_stims), function(i) afni_command_switch(afni_stims[[i]], i, "iresp")))
  
  #browser()
  
  if ( (length(opt_stim_times) + length(opt_stim_times_IM)) > 0) {
    ## if we use `afni_hrf` that use -stim_times, then we use local times
    global_times <- FALSE
  } else {
    ## otherwise global_times is irrelevant, since values rather than times are provided.
    global_times <- TRUE
  }
  
  ## support pre-refactor datasets that stored scans/mask directly
  ## When using nodata, we don't need actual input files
  inputs <- if (!is.null(opts[["nodata"]])) {
              NULL  # No input files needed for nodata mode
            } else if (!is.null(dataset$scans)) {
              dataset$scans
            } else if (!is.null(dataset$backend$source)) {
              dataset$backend$source
            } else {
              stop("No input scans found in dataset")
            }

  mask_file <- if (!is.null(opts[["nodata"]])) {
                 NULL  # No mask needed for nodata mode
               } else if (!is.null(dataset$mask_file)) {
                 dataset$mask_file
               } else if (!is.null(dataset$backend$mask_source)) {
                 dataset$backend$mask_source
               } else {
                 NULL
               }

  # For multiple runs with nodata, create concat string
  concat_str <- NULL
  # Get blockids from the event model
  event_blockids <- model$event_model$blockids
  
  if (!is.null(event_blockids) && length(unique(event_blockids)) > 1 && !is.null(opts[["nodata"]])) {
    # Calculate run start points based on nodata parameter
    # nodata[1] is total volumes, nodata[2] is TR
    # Need to calculate volumes per run
    n_runs <- length(unique(event_blockids))
    volumes_per_run <- opts[["nodata"]][1] / n_runs
    
    # Generate run starts: 0, volumes_per_run, 2*volumes_per_run, ...
    run_starts <- seq(0, by = volumes_per_run, length.out = n_runs)
    concat_str <- paste0("1D: ", paste(as.integer(run_starts), collapse=" "))
  }
  
  cmdlines <- list(nodata=if (!is.null(opts[["nodata"]])) paste(opts[["nodata"]], collapse=" ") else NULL,
                   concat=concat_str,
                   input=if (!is.null(inputs)) paste0(inputs) else NULL,
                   mask=if (!is.null(mask_file)) paste0(mask_file) else NULL,
                   polort=if (opts[["polort"]] > 0) opts[["polort"]] else -1,
                   global_times=if (global_times) TRUE else NULL,
                   num_stimts=length(afni_stims),
                   num_glt=length(gltfiles),
                   stim_file=opt_stim_files,
                   stim_label=opt_stim_labels,
                   ortvec=opt_stim_ortvecs,
                   censor=if (!is.null(opts[["censor"]])) "censor.1D" else NULL,
                   stim_times=opt_stim_times,
                   stim_times_IM=opt_stim_times_IM,
                   TR_times=opts[["TR_times"]],
                   iresp=opt_stim_iresp,
                   gltsym=lapply(seq_along(gltfiles), function(i) paste(gltfiles[i], collapse=" ")),
                   glt_label=lapply(seq_along(gltnames), function(i) paste(i, gltnames[i], collapse=" ")),
                   nofullf_atall=opts[["nofullf_atall"]],
                   fout=opts[["fout"]],
                   rout=opts[["rout"]],
                   tout=opts[["tout"]],
                   bout=opts[["bout"]],
                   noFDR=opts[["noFDR"]],
                   cbucket=opts[["cbucket"]],
                   bucket=opts[["bucket"]],
                   nocond=opts[["nocond"]],
                   x1D_stop=opts[["x1D_stop"]],
                   jobs=opts[["jobs"]],
                   errts=if (!is.null(opts[["errts"]])) opts[["errts"]] else NULL,
                   float=TRUE)
  
  cmd <- .make_decon_command_str(cmdlines)
  
  list(cmd=cmd, cmdlines=cmdlines, afni_stims=afni_stims, afni_baseline_mats=afni_baseline_mats,
       gltfiles=gltfiles, gltnames=gltnames, glts=glts, gltstr=gltstr, censor=opts$censor)
}


#' @importFrom assertthat assert_that
NULL