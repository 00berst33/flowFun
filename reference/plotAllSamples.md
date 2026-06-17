# plotAllSamples

Plot chosen gate for all samples

## Usage

``` r
plotAllSamples(gs, xdim, ydim, subset, node)
```

## Arguments

- gs:

  A `GatingSet`

- xdim:

  The channel to plot on the x-axis

- ydim:

  The channel to plot on the y-axis

- subset:

  The population whose cells should be displayed on the plot

- node:

  The gate to be visualized

  Note that to use this function, `gs` must have at least one gate.

## Value

A `ggplot` drawn with `ggcyto`
