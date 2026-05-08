# getSampleMetaclusterMFIs

Get a table of MFIs, where each row is a sample and each column is a
metacluster.

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

A table, where rows are samples and columns are metaclusters.
