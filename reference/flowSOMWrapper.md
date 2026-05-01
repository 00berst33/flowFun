# flowSOMWrapper

flowSOMWrapper

## Usage

``` r
flowSOMWrapper(
  input,
  cols_to_cluster,
  num_clus,
  seed = NULL,
  fsom_file = NULL,
  ...
)
```

## Arguments

- input:

  A `flowSet`

- cols_to_cluster:

  A numeric or character vector specifying which columns should be used
  for clustering.

- num_clus:

  The number of metaclusters to create.

- seed:

  Optional, a seed for reproducibility. Default is `NULL`.

- fsom_file:

  An .rds filename. If not `NULL` (default), the resulting FlowSOM
  object will be saved under this name.

- ...:

  Additional parameters to pass to
  [`FlowSOM()`](https://rdrr.io/pkg/FlowSOM/man/FlowSOM.html).

## Value

A data.table annotated with cluster and metacluster assignments.
