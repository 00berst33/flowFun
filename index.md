# flowFun

R package integrating an end-to-end pipeline for the analysis of high
parameter cytometry data. Functions for pre-processing, cell type
identification, differential analysis, data visualization, and more are
included.

## Getting Started

Vignettes and other helpful documentation may be viewed on the package’s
main page: <https://00berst33.github.io/flowFun>

## Installation

This package can be installed with vignettes from the RStudio console as
follows:

``` r

install.packages("BiocManager")
install.packages("devtools")

# Install with vignettes
#   The package 'flowFunData' must first be installed
devtools::install_github("00berst33/flowFunData")
#   Then install flowFun
devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = TRUE)
```

Or to install without vignettes:

``` r

install.packages("devtools")

# Install without vignettes
#   No installation of 'flowFunData', set argument `build_vignettes` to FALSE
devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = FALSE)
```
