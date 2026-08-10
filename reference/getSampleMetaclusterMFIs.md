# getSampleMetaclusterMFIs

Generate a table of MFIs from clustered data, where each row is a sample
and each column is a metacluster.

## Usage

``` r
getSampleMetaclusterMFIs(
  input,
  col,
  populations = NULL,
  transformation = NULL,
  inverse = FALSE
)
```

## Arguments

- input:

  A `data.table` or `GatingSet`

- col:

  The channel to find sample-metacluster MFIs for.

- populations:

  A `character` list of the names of metaclusters/populations to
  calculate MFIs for.

- transformation:

  A
  [`transformList`](https://rdrr.io/pkg/flowCore/man/transformList-class.html)
  specifying the transformation applied to the data, if any. Required
  when `input` is a `data.table` and `inverse = TRUE`

- inverse:

  Boolean, whether data should be back-transformed before calculating
  MFIs. Only valid when `input` is a `GatingSet` containing a
  transformation.

## Value

A `data.table`, where rows are samples and columns are metaclusters.

## Examples

``` r
## data.table method
  file <- system.file("extdata", "fsom_dt.rds", package = "flowFunData")
  table <- readRDS(file)

  # define populations of interest
  pop_of_interest <- c("CD8 T cells", "CD4 T cells", "NK T cells", "gdT cells", "B cells")

  mfis <- getSampleMetaclusterMFIs(table,
                                   col = "FITC-A",
                                   populations = pop_of_interest)
#> Error in getSampleMetaclusterMFIs.data.table(table, col = "FITC-A", populations = pop_of_interest): argument "transformation" is missing, with no default

## GatingSet method
  # not run: assume gs is a GatingSet

  # get all children nodes/subpopulations of "live" population
  subpops <- gs_pop_get_children(gs, "live")
#> Error in gs_pop_get_children(gs, "live"): could not find function "gs_pop_get_children"

  mfis <- getSampleMetaclusterMFIs(gs,
                                   col = "FITC-A",
                                   populations = subpops,
                                   inverse = TRUE) # no need to specify `transformation` parameter
#> Error: object 'gs' not found

## Save table
write.csv(mfis, "mfi_table_fitc.csv")
#> Error in eval(expr, p): object 'mfis' not found
```
