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

## Usage

```r
library(afnireg)

# Create an AFNI HRF specification
hrf_spec <- afni_hrf(type = "SPMG1")

# Build a model for AFNI
model <- afni_lm(
  dataset = my_data,
  event_model = my_events,
  hrf_spec = hrf_spec
)

# Run the model in AFNI
results <- run(model)
```

## Documentation

See the package vignette for a complete tutorial:
```r
vignette("afni_tutorial", package = "afnireg")
```

## License

GPL (>= 2)
