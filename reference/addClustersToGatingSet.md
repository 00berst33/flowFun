# addClustersToGatingSet

Uses cluster labels to add corresponding gates to the GatingSet

## Usage

``` r
addClustersToGatingSet(table, gs, parent_gate, fsom_file = NULL)
```

## Arguments

- table:

  A data.table or tibble containing expression data for all samples.
  Must have a column named `.id` indicating sample ID, and a column
  `Metacluster` with cluster labels.

- gs:

  The `GatingSet` to add gates to.

- parent_gate:

  `character` indicating which population is the parent of the clusters.

- fsom_file:

  Optional, a character giving a path to an .rds file for a FlowSOM
  object. Default is `NULL`, which checks for a filename in
  `attr(table, "fsom_filename")`

  This function should be used after a satisfactory clustering has been
  obtained using the human-in-the-loop approach outlined by the
  vignette. The resulting table should contain expression data for each
  sample, with each cell assigned a cluster label. These labels will
  then be used to construct a gate for each cluster.
