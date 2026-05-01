# plotUMAP

Plots UMAP colored by metacluster.

## Usage

``` r
plotUMAP(input, num_cells = 5000, labels = NULL, colors = NULL, seed = NULL)
```

## Arguments

- input:

  A FlowSOM object or data table.

- num_cells:

  The number of cells to use for the dimension reduction. Default is
  5000.

- labels:

  Optional, a vector specifying the order in which metacluster name will
  appear in the legend. All metacluster names present in `input` should
  exists in this vector.

- colors:

  Optional, a vector of colors to use in the plot for each metacluster.

- seed:

  Optional, a seed for reproducibility.

## Value

A UMAP plot drawn with ggplot2.

## Examples

``` r
file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
fsom <- FlowSOM::FlowSOM(file,
                         colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
                         nClus = 10,
                         seed = 42,
                         xdim = 6,
                         ydim = 6)
#> Error in if (file.info(input[i])$isdir) {    toAdd <- c(toAdd, list.files(input[i], pattern = pattern,         full.names = TRUE))    toRemove <- c(toRemove, i)}: missing value where TRUE/FALSE needed
#> Timing stopped at: 0.001 0 0.001

plotUMAP(fsom)
#> Error: object 'fsom' not found

plotUMAP(fsom)
#> Error: object 'fsom' not found
```
