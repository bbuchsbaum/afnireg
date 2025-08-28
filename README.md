# afnireg

AFNI Integration for fMRI Analysis in R

## Overview

The `afnireg` package provides an R interface to AFNI's 3dDeconvolve program for fMRI analysis. It allows users to specify complex experimental designs and HRF models in R using the fmridesign framework and execute them using AFNI's optimized estimation routines.

## Installation

```r
# Install from GitHub
remotes::install_github("bbuchsbaum/afnireg")
```

## Dependencies

This package depends on:
- fmridesign: Event model specification
- fmrihrf: HRF modeling
- fmridataset: Dataset handling

## Features

- Native AFNI HRF specifications
- Trialwise modulation support
- Automatic contrast conversion to GLT format
- Support for multiple HRF basis functions
- Multi-run designs
- Censoring and nuisance regression

## Quick Start

For a comprehensive tutorial on using `afnireg` to translate fMRI models to AFNI's 3dDeconvolve format, see our [Getting Started vignette](https://bbuchsbaum.github.io/afnireg/articles/afni_tutorial.html).

## Usage

```r
library(afnireg)
library(fmridesign)
library(fmrihrf)

# 1) Specify an event and baseline model in R
TR <- 2
sframe <- sampling_frame(blocklens = 140, TR = TR)

emodel <- event_model(
  onset ~ hrf(stim),
  data = events_df,       # data.frame with columns: onset, stim, run, ...
  block = ~ run,
  sampling_frame = sframe
)

bmodel <- baseline_model(basis = "bs", degree = 5, sframe = sframe)
fmodel <- fmri_model(emodel, bmodel)

# 2a) With scans: build an AFNI spec
dset <- fmridataset::fmri_dataset(
  scans = c("run1.nii.gz", "run2.nii.gz"),
  mask  = "mask.nii.gz",
  TR    = TR,
  run_length = 140,
  event_table = events_df
)

alm <- afni_lm(fmodel, dset, options = list(bucket = "stats_afni"))
run(alm, outdir = "glm_afni_output")

# 2b) Or do a dry run without data (writes .xmat.1D and script)
alm_nodata <- afni_lm(fmodel, dataset = NULL, nodata = c(140, TR), x1D_stop = TRUE)
run(alm_nodata, outdir = "glm_afni_x1d", execute = FALSE)

# AFNI-native HRF example (uses -stim_times)
emodel_afni <- event_model(
  onset ~ afni_hrf(stim, basis = "block", durations = 2),
  data = events_df, block = ~ run, sampling_frame = sframe
)
fmodel_afni <- fmri_model(emodel_afni, bmodel)
alm_afni <- afni_lm(fmodel_afni, dset)
```

## Documentation

### Online Documentation
- [Package documentation](https://bbuchsbaum.github.io/afnireg/)
- [Getting Started: Translating to AFNI 3dDeconvolve](https://bbuchsbaum.github.io/afnireg/articles/afni_tutorial.html)

### Vignettes
View the vignette locally after installation:
```r
vignette("afni_tutorial", package = "afnireg")
```

## License

GPL (>= 2)
