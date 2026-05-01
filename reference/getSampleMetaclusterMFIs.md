# getSampleMetaclusterMFIs

Get a table of MFIs, where each row is a sample and each column is a
metacluster.

## Usage

``` r
getSampleMetaclusterMFIs(input, col, sample_df, meta_to_use = NULL)
```

## Arguments

- input:

  description

- col:

  The channel to find sample-metacluster MFIs for.

- sample_df:

  If `input` is a FlowSOM object, a data frame from
  [`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md).

- meta_to_use:

  A `character` list of the names of metaclusters/populations to
  calculate MFIs for. Default is all.

## Value

A table, where rows are samples and columns are metaclusters.
