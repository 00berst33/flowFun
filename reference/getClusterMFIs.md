# getClusterMFIs

getClusterMFIs

## Usage

``` r
getClusterMFIs(input, cols_to_use = NULL)
```

## Arguments

- input:

  A data frame or table, with column `"Cluster"`.

- cols_to_use:

  A vector specifying which columns to calculate medians for. Default is
  all columns used for clustering.

## Value

A data frame, where rows are clusters and columns are channels.
