#' AFNI Linear Model Functions
#'
#' This file contains the main interface functions for creating and running
#' AFNI linear models using 3dDeconvolve.
#'

#' Generate AFNI command for fmri_config
#' @keywords internal
#' @noRd
gen_afni_lm.fmri_config <- function(x, ...) {
  mspec <- gen_model_spec(x)
  gen_afni_lm.fmri_model(mspec,...)
  #gen_afni_lm(x$event_model, x$baseline_model, ...)
}

#' Generate AFNI command for fmri_model
#' @keywords internal
#' @noRd
#' @importFrom fmridesign baseline_model event_model
#' @importFrom fmrireg fmri_model
gen_afni_lm.fmri_model <- function(x, dataset, working_dir=".", 
                                    polort=-1, 
                                    jobs=1, 
                                    censor=NULL, 
                                    iresp=FALSE, 
                                    nodata=NULL,
                                    TR_times=x$sampling_frame$TR,
                                    x1D_stop=FALSE,
                                    nofullf_atall=TRUE,
                                    ...) {
  
  dset <- if (!is.null(dataset)) {
    dataset
  } else {
    # Handle case where dataset might be embedded in the model
    if (!is.null(x$dataset)) {
      x$dataset
    } else {
      NULL
    }
  }
  
  if (is.null(dset) && is.null(nodata)) {
    stop("Either 'dataset' or 'nodata' must be provided")
  }
  
  options <- list(...)
  # Set defaults for common options, but preserve any existing values in options
  if (is.null(options$polort)) options$polort <- polort
  if (is.null(options$jobs)) options$jobs <- jobs
  if (is.null(options$censor)) options$censor <- censor
  if (is.null(options$iresp)) options$iresp <- iresp
  if (is.null(options$nodata)) options$nodata <- nodata
  if (is.null(options$TR_times)) options$TR_times <- TR_times
  if (is.null(options$x1D_stop)) options$x1D_stop <- x1D_stop
  if (is.null(options$nofullf_atall)) options$nofullf_atall <- nofullf_atall
  
  cmd <- build_decon_command(x, dset, working_dir, options)
  
  ret <- list(
    model = x,
    dataset = dset,
    working_dir = working_dir,
    options = options,
    cmd = cmd
  )
  
  class(ret) <- "afni_lm_spec"
  ret
}


#' create an afni_lm spec
#' 
#' @param fmri_mod the fmri_model
#' @param dataset the dataset
#' @param working_dir the working directory
#' @param polort polynomial detrending order
#' @param jobs number of parallel jobs
#' @param censor censoring vector
#' @param iresp whether to output IRFs
#' @param nodata nodata option for testing
#' @param TR_times TR in seconds
#' @param x1D_stop whether to stop after creating design matrix
#' @param nofullf_atall skip full F-test
#' @param options optional list of additional options that override defaults
#' @param ... additional 3dDeconvolve options
#' 
#' @export
afni_lm <- function(fmri_mod, dataset, working_dir=".", polort=-1, jobs=1, 
                   censor=NULL, iresp=FALSE, nodata=NULL,
                   TR_times=fmri_mod$sampling_frame$TR,
                   x1D_stop=FALSE,
                   nofullf_atall=TRUE,
                   options=NULL,
                   ...) {
  
  # If options list is provided, use those values to override defaults
  if (!is.null(options) && is.list(options)) {
    if (!is.null(options$polort)) polort <- options$polort
    if (!is.null(options$jobs)) jobs <- options$jobs
    if (!is.null(options$censor)) censor <- options$censor
    if (!is.null(options$iresp)) iresp <- options$iresp
    if (!is.null(options$nodata)) nodata <- options$nodata
    if (!is.null(options$TR_times)) TR_times <- options$TR_times
    if (!is.null(options$x1D_stop)) x1D_stop <- options$x1D_stop
    if (!is.null(options$nofullf_atall)) nofullf_atall <- options$nofullf_atall
  }
  
  gen_afni_lm.fmri_model(fmri_mod, dataset, working_dir, polort, jobs, 
                         censor, iresp, nodata, TR_times, x1D_stop, 
                         nofullf_atall, ...)
}

#' Construct an afni_lm_spec object
#'
#' Creates an `afni_lm_spec` object bundling the model, dataset, options,
#' working directory, and generated 3dDeconvolve command.
#'
#' @param model An `fmri_model` used to generate the AFNI command
#' @param dataset An `fmri_dataset` or NULL if using `nodata`
#' @param working_dir Character path used when generating files
#' @param options A named list of options passed to command generation
#' @param cmd A list with `cmd` (string) and other command metadata
#' @return An object of class `afni_lm_spec`
#' @export
afni_lm_spec <- function(model, dataset, working_dir, options, cmd) {
  structure(
    list(model = model, dataset = dataset, working_dir = working_dir, 
         options = options, cmd = cmd),
    class = "afni_lm_spec"
  )
}

#' Print method for afni_lm_spec
#' 
#' @param x An afni_lm_spec object
#' @param ... Additional arguments
#' @export
print.afni_lm_spec <- function(x,...) {
  cat("afni_lm_spec", "\n")
  cat("  model type: fmri_model", "\n")
  cat("  working directory: ", x$working_dir, "\n")
  cat("  dataset: ", if(!is.null(x$dataset)) "provided" else "none", "\n")
  cat("  options: ", "\n")
  print(x$options)
}

#' Generate AFNI command for afni_lm_spec
#' @keywords internal
#' @noRd
gen_afni_lm.afni_lm_spec <- function(x, ...) {
  x
}

#' Generic run function
#' 
#' @param x An object to run
#' @param ... Additional arguments
#' @export
run <- function(x, ...) {
  UseMethod("run")
}

#' Run AFNI Linear Model Specification
#' 
#' Execute an AFNI linear model specification using 3dDeconvolve.
#' 
#' @param x An afni_lm_spec object
#' @param outdir Output directory
#' @param execute Logical; execute the command
#' @param execfun Function to execute the command (default: system)
#' @param prepend String to prepend to command
#' @param ... Additional arguments
#' @return Result of execution
#' @export
run.afni_lm_spec <- function(x, outdir, execute=TRUE, execfun=system, prepend="",...) {
  start_dir <- getwd()
  res <- try({
    if (!file.exists(outdir)) {
      dir.create(outdir)
    } else {
      warning(paste("glm output directory: ", outdir, " already exists"))
      outdir <- next_dir_name(outdir)
      dir.create(outdir)
      warning(paste("outputting to: ", outdir))
    }
    print(paste("setting directory:", outdir))
    setwd(outdir)
    
    write_stim_files(x$cmd$afni_stims)
    
    if (!is.null(x$cmd$gltstr)) {
      write_glts(x$cmd$gltstr, x$cmd$gltfiles)
    }
    
    if (!is.null(x$cmd$afni_baseline_mats)) {
      write_baseline_mats(x$cmd$afni_baseline_mats)
    }
    
    if (!is.null(x$cmd$censor)) {
      write_censor_file(".", x$cmd$censor)
    }
    
    #if (reml) {
    #  x$cmd$cmd <- paste(x$cmd$cmd, "-x1D_stop")
    #}
    
    write(x$cmd$cmd, "3ddeconvolve.sh")
    
    if (execute) {
      execfun(paste(prepend, x$cmd$cmd))
      
      #if (reml) {
      #  execfun(paste0("./", x$options$bucket, ".REML_cmd"))
      #}
    }
  })
  
  setwd(start_dir)
}
