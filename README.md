
# shinySTRegister

<!-- badges: start -->
<!-- badges: end -->

The goal of shinySTRegister is to provide a simple (if barebones) interface for
registering (i.e. aligning) two spatial transcriptomics samples. This is done by
specifying three points on the first sample (to be aligned) and their
corresponding projections on the second (reference) sample.

The corresponding affine transformation is then determined and applied to the
first sample, which can be saved as RDS for future use.

## Installation

You can install the development version of shinySTRegister from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JonSulc/shinySTRegister")
```

## Example

Given two samples (to_align and reference) loaded as SpatialExperiment objects,
the interface can be started with:

``` r
library(shinySTRegister)
shiny_register(to_align, reference, nbins = 400)
```

The `nbins` parameter controls the resolution of the image, counts being binned
into a grid of $nbins \times nbins$ pixels. This only affects the display
resolution, the transformation being applied to the raw data coordinates.

## Convenience functions

To help with the setup for this, there are a few convenience functions for the
loading and display of data from certain technologies.
`load_as_spatial_experiment` takes a counts filename (spot/cell $\times$
feature) and a spot/cell metadata filename or data.frame to create the
`SpatialExperiment` object. The `id_columns` argument specifies the columns
identifying the corresponding spots or cells to combine the data. If `NULL`,
they are assumed to be in the 
