# getMetaclusterMFIs

getMetaclusterMFIs

## Usage

``` r
getMetaclusterMFIs(input, cols_to_use = NULL)
```

## Arguments

- input:

  A data frame or table, with column `"Metacluster"`.

- cols_to_use:

  A vector specifying which columns to calculate medians for. Default is
  all columns used for clustering.

## Value

A data frame, where rows are metaclusters and columns are channels.
