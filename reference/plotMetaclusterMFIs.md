# plotMetaclusterMFIs

Plot a heatmap of metacluster MFIs.

## Usage

``` r
plotMetaclusterMFIs(input, cols_to_use = NULL, ...)
```

## Arguments

- input:

  A FlowSOM object or data table.

- cols_to_use:

  A character vector specifying which markers you would like to include
  in the plot. Default is the markers used for clustering.

- ...:

  Additional parameters to pass to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

## Value

A heatmap of MFIs made with
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
where each column is a marker of interest, and each row is a
metacluster.

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
#> Timing stopped at: 0.002 0 0.002

plotMetaclusterMFIs(fsom)
#> Error: object 'fsom' not found
```
