pkgname <- "afnireg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('afnireg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("afni_hrf")
### * afni_hrf

flush(stderr()); flush(stdout())

### Name: afni_hrf
### Title: construct an native AFNI hrf specification for '3dDeconvolve'
###   with the 'stim_times' argument.
### Aliases: afni_hrf

### ** Examples

# Create SPM canonical HRF specification
hrf1 <- afni_hrf(onsets = c(10, 30, 50), basis = "spmg1")

# Create block HRF with duration
hrf2 <- afni_hrf(onsets = c(10, 30, 50), durations = 5, basis = "block")

# Create tent basis HRF
hrf3 <- afni_hrf(onsets = c(10, 30, 50), basis = "tent", 
                 start = 0, stop = 20, nbasis = 10)



cleanEx()
nameEx("afni_trialwise")
### * afni_trialwise

flush(stderr()); flush(stdout())

### Name: afni_trialwise
### Title: construct a native AFNI hrf specification for '3dDeconvolve' and
###   individually modulated events using the 'stim_times_IM' argument.
### Aliases: afni_trialwise

### ** Examples


tw <- afni_trialwise("trialwise", basis="gamma", onsets=seq(1,100,by=5))




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
