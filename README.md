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
