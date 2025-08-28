#' AFNI HRF Specifications
#'
#' This file contains all AFNI HRF (Hemodynamic Response Function) specifications,
#' constructors, and related functions for creating AFNI-compatible HRF models.
#'

#' AFNI HRF Constructor Function
#'
#' @description
#' The `AFNI_HRF` function creates an object representing an AFNI-specific hemodynamic response function (HRF). It is a class constructor for AFNI HRFs.
#'
#' @param name A string specifying the name of the AFNI HRF.
#' @param nbasis An integer representing the number of basis functions for the AFNI HRF.
#' @param params A list containing the parameter values for the AFNI HRF.
#' @param span A numeric value representing the span in seconds of the HRF. Default is 24.
#'
#' @return An AFNI_HRF object with the specified properties.
#'
#' @seealso HRF
#'
#' @export
#' @rdname AFNI_HRF-class
AFNI_HRF <- function(name, nbasis, params, span = 24) {
  structure(name,
            nbasis = as.integer(nbasis),
            params = params,
            span = span,
            class = c("AFNI_HRF", "HRF"))
  
}


#' @export
as.character.AFNI_HRF <- function(x,...) {
  params <- attr(x, "params")
  param_str <- if (length(params) > 0) {
    paste(unlist(params), collapse=",")
  } else {
    ""
  }
  paste0(x, "(", param_str, ")")
}

#' construct an native AFNI hrf specification for '3dDeconvolve' with the 'stim_times' argument.
#' 
#' @param ... Variables to include in the HRF specification
#' @param basis Character string specifying the basis function type
#' @param onsets Numeric vector of event onset times
#' @param durations Numeric vector of event durations
#' @param prefix Character string prefix for the term
#' @param subset Expression for subsetting events
#' @param nbasis Number of basis functions
#' @param contrasts Contrast specifications
#' @param id Character string identifier for the term
#' @param lag Numeric lag in seconds
#' @param precision Numeric precision for convolution
#' @param summate Logical whether to summate overlapping responses
#' @param start the start of the window for sin/poly/csplin models
#' @param stop the stop time for sin/poly/csplin models
#' @return an \code{afni_hrfspec} instance with class "afni_hrfspec"
#' @examples
#' # Create SPM canonical HRF specification
#' hrf1 <- afni_hrf(onsets = c(10, 30, 50), basis = "spmg1")
#' 
#' # Create block HRF with duration
#' hrf2 <- afni_hrf(onsets = c(10, 30, 50), durations = 5, basis = "block")
#' 
#' # Create tent basis HRF
#' hrf3 <- afni_hrf(onsets = c(10, 30, 50), basis = "tent", 
#'                  start = 0, stop = 20, nbasis = 10)
#' @export
afni_hrf <- function(..., basis=c("spmg1", "block", "dmblock",
                                  "tent",   "csplin", "poly",  "sin", "sine",
                                  "gam", "gamma", "spmg2", "spmg3", "wav"),
                                  onsets=NULL, durations=NULL, prefix=NULL, subset=NULL,
                                  nbasis=1, contrasts=NULL, id=NULL,
                                  lag=0, precision = 0.3, summate = TRUE,
                                  start=NULL, stop=NULL) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- tolower(basis)
  if (basis == "sin") basis <- "sine"
  if (basis == "gam") basis <- "gamma"
  basis <- match.arg(basis)
  
  vars <- as.list(base::substitute(list(...)))[-1]
  # Convert raw expressions to quosures for compatibility with construct_event_term
  vars_quos <- lapply(vars, function(expr) rlang::new_quosure(expr, env = rlang::caller_env()))
  # Requires parse_term (assuming it exists elsewhere now or is removed)
  # parsed <- parse_term(vars, "afni_hrf") 
  # term <- parsed$term
  # label <- parsed$label
  # --- Need to replace parse_term dependency --- 
  # Simplified naming based on input expressions/symbols for now
  var_labels <- sapply(vars, rlang::as_label)
  varnames <- sapply(vars, function(v) {
       if (rlang::is_symbol(v)) as.character(v) else make.names(rlang::as_label(v))
  })
  term <- vars_quos # Store quosures instead of raw expressions
  label <- paste0("afni_hrf(", paste0(var_labels, collapse=","), ")")
  # --- End replacement --- 
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_hrf does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=nbasis, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=nbasis, b=start, c=stop)
  }
  
  
  # varnames <- if (!is.null(prefix)) {
  #   paste0(prefix, "_", term)
  # } else {
  #   term
  # } 
  # Use the new varnames logic
  if (!is.null(prefix)) {
      varnames <- paste0(prefix, "_", varnames)
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  
  
  # Handle contrasts parameter
  cset <- if (!is.null(contrasts) && inherits(contrasts, "contrast_spec")) {
    # Single contrast - wrap in a list
    list(contrasts)
  } else if (!is.null(contrasts) && inherits(contrasts, "contrast_set")) {
    # Already a contrast set - it's a list with contrast_set class
    # The contrast_set itself IS the list of contrasts
    contrasts
  } else if (!is.null(contrasts) && is.list(contrasts)) {
    # List of contrasts (check if they're all contrast_spec objects)
    if (all(sapply(contrasts, inherits, "contrast_spec"))) {
      contrasts
    } else {
      NULL
    }
  } else { 
    NULL 
  } # Default to NULL if not correct type
  
  ret <- list(
    name=termname,
    id=id,
    varnames=varnames,
    vars=term, # Store original expressions/symbols
    label=label,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=rlang::enexpr(subset), # Capture subset expression
    lag=lag,
    precision = precision,
    summate = summate,
    contrasts=cset)
  
  class(ret) <- c("afni_hrfspec", "hrfspec", "list")
  ret
  
}

#' construct a native AFNI hrf specification for '3dDeconvolve' and individually modulated events using the 'stim_times_IM' argument.
#' 
#' @param label name of regressor
#' @param basis Character string specifying the basis function type
#' @param onsets Numeric vector of event onset times
#' @param durations Numeric vector of event durations (default 0)
#' @param subset Expression for subsetting events
#' @param id Character string identifier for the term
#' @param precision Numeric precision for convolution (default 0.3)
#' @param summate Logical whether to summate overlapping responses (default TRUE)
#' @param start start of hrf (for multiple basis hrfs)
#' @param stop end of hrf (for multiple basis hrfs)
#' 
#' @examples
#' 
#' tw <- afni_trialwise("trialwise", basis="gamma", onsets=seq(1,100,by=5))
#' 
#' @export
#' @return an \code{afni_trialwise_hrfspec} instance
afni_trialwise <- function(label, basis=c("spmg1", "block", "dmblock", "gamma", "wav"),
                     onsets=NULL, durations=0, subset=NULL,
                      id=NULL, start=0, stop=22,
                      precision = 0.3, summate = TRUE) {
  
  # Capture label as a quosure to prevent evaluation
  label_quo <- rlang::enquo(label)
  
  # Extract the variable name as a string
  label_name <- rlang::as_label(label_quo)
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  basis <- tolower(basis)
  if (basis == "gam") basis <- "gamma"
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_trialwise does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=1, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=1, b=start, c=stop)
  }
  
  
  if (is.null(id)) {
    id <- label_name
  }  
  
  ret <- list(
    name=label_name,
    varname=label_name,
    id=id,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    subset=rlang::enexpr(subset),
    precision = precision,
    summate = summate)
  
  class(ret) <- c("afni_trialwise_hrfspec", "hrfspec", "list")
  ret
  
}


## AFNI HRF Constructor Functions
#'
#' AFNI-native HRF basis constructors. These functions create `AFNI_HRF`
#' objects corresponding to AFNI’s built-in basis families.
#'
#' @param d Duration parameter (seconds) used by SPMG/BLOCK/WAV families
#' @param p Penalty or additional parameter for BLOCK/dmBLOCK
#' @param b Start time (seconds) for multi-basis families (e.g., TENT/CSPLIN/POLY/SIN)
#' @param c Stop time (seconds) for multi-basis families (e.g., TENT/CSPLIN/POLY/SIN)
#' @param n Number of basis functions/knots/degrees for multi-basis families
#' @param q Shape parameter for `AFNI_GAM`
#'
#' @section Constructors:
#' - `AFNI_SPMG1(d=1)` – SPMG1 (canonical), duration `d` seconds
#' - `AFNI_SPMG2(d=1)` – SPMG2, duration `d`
#' - `AFNI_SPMG3(d=1)` – SPMG3, duration `d`
#' - `AFNI_BLOCK(d=1, p=1)` – BLOCK, duration `d`, penalty `p`
#' - `AFNI_dmBLOCK(d=1, p=1)` – dmBLOCK, duration `d`, penalty `p`
#' - `AFNI_TENT(b=0, c=18, n=10)` – TENT from `b` to `c` with `n` knots
#' - `AFNI_CSPLIN(b=0, c=18, n=6)` – CSPLIN from `b` to `c` with `n` splines
#' - `AFNI_POLY(b=0, c=18, n=10)` – POLY from `b` to `c` with `n` degrees
#' - `AFNI_SIN(b=0, c=18, n=10)` – SIN from `b` to `c` with `n` harmonics
#' - `AFNI_GAM(p=8.6, q=.547)` – GAM with AFNI default parameters
#' - `AFNI_WAV(d=1)` – WAV (duration parameter)
#'
#' @name AFNI_HRF_constructors
NULL

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_SPMG1 <- function(d=1) AFNI_HRF(name="SPMG1", nbasis=as.integer(1), params=list(d=d)) 

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_SPMG2 <- function(d=1) AFNI_HRF(name="SPMG2", nbasis=as.integer(2), params=list(d=d))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_SPMG3 <- function(d=1) AFNI_HRF(name="SPMG3", nbasis=as.integer(3), params=list(d=d))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_BLOCK <- function(d=1,p=1) AFNI_HRF(name="BLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_dmBLOCK <- function(d=1,p=1) AFNI_HRF(name="dmBLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_TENT <- function(b=0,c=18, n=10) AFNI_HRF(name="TENT", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_CSPLIN <- function(b=0,c=18, n=6) AFNI_HRF(name="CSPLIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_POLY <- function(b=0,c=18, n=10) AFNI_HRF(name="POLY", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_SIN <- function(b=0,c=18, n=10) AFNI_HRF(name="SIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_GAM <- function(p=8.6,q=.547) AFNI_HRF(name="GAM", nbasis=as.integer(1), params=list(p=p,q=q))

#' @rdname AFNI_HRF_constructors
#' @export
AFNI_WAV <- function(d=1) AFNI_HRF(name="WAV", nbasis=as.integer(1), params=list(d=1))


#' @keywords internal
#' @noRd
get_AFNI_HRF <- function(name, nbasis=1, duration=1, b=0, c=18) {
  name <- toupper(name)
  if (name == "SINE") name <- "SIN"
  if (name == "GAMMA") name <- "GAM"
  
  switch(name,
         "SPMG1" = AFNI_SPMG1(d=duration),
         "SPMG2" = AFNI_SPMG2(d=duration),
         "SPMG3" = AFNI_SPMG3(d=duration),
         "BLOCK" = AFNI_BLOCK(d=duration, p=1),
         "DMBLOCK" = AFNI_dmBLOCK(d=duration, p=1),
         "TENT" = AFNI_TENT(b=b, c=c, n=nbasis),
         "CSPLIN" = AFNI_CSPLIN(b=b, c=c, n=nbasis),
         "POLY" = AFNI_POLY(b=b, c=c, n=nbasis),
         "SIN" = AFNI_SIN(b=b, c=c, n=nbasis),
         "GAM" = AFNI_GAM(),
         "WAV" = AFNI_WAV(d=1),
         stop(paste("unsupported HRF type: ", name)))
}


#' @importFrom assertthat assert_that
#' @importFrom rlang enexpr enquo as_label
NULL
