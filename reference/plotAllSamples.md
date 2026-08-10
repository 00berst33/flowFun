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

## Value

A `ggplot` drawn with `ggcyto`

## Details

Intended for checking the results of newly added gates, typically after
a call to
[`gt_gating()`](https://rdrr.io/pkg/openCyto/man/gt_gating.html). Note
that to use this function, `gs` must have at least one gate.

## Examples

``` r
# load in data
files <- system.file("extdata", "samples", package = "flowFunData")
cs <- flowWorkspace::load_cytoset_from_fcs(path=dir, pattern=".fcs", truncate_max_range = TRUE)
#> Error in dir(path, pattern, full.names = TRUE): invalid 'path' argument

# Make GatingSet
gs <- flowWorkspace::GatingSet(cs)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'GatingSet': object 'cs' not found

# Generate default gating scheme
gt_table <- generateGatingTable(gs, collapse_data = TRUE)
#> Error: object 'gs' not found
gt <- openCyto::gatingTemplate(gt_table)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'gatingTemplate': object 'gt_table' not found

# Apply gating scheme
  # not run: assume transformations and compensations were performed as needed
openCyto::gt_gating(gt, gs)
#> Error: object 'gt' not found

## Plotting results
# View nonDebris gate
plotAllSamples(gs, "FSC-A", "SSC-A", "nonMargins", "nonDebris")
#> Error: object 'gs' not found

# View singlets gate
plotAllSamples(gs, "FSC-A", "FSC-H", "nonDebris", "singlets")
#> Error: object 'gs' not found
```
