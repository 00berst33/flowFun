# plotMarkerDensityByChannel

plotMarkerDensityByChannel

## Usage

``` r
plot1DMarkerDensities(
  gs,
  channel,
  population = "root",
  facet_by = c("samples", "subpopulations"),
  inverse = FALSE,
  verbose = TRUE
)
```

## Arguments

- gs:

  A `GatingSet`

- channel:

  The channel to plot densities for

- population:

  A `character` indicating which subpopulation to use

- inverse:

  A boolean, whether or not the data should be inverse transformed
  before plotting. `FALSE` by default

- verbose:

  Boolean specifying whether or not to print progress updates as
  function runs. Default is `TRUE`

## Value

A faceted plot with density plots of the chosen channel for each sample
