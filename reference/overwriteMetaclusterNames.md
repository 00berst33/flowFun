# overwriteMetaclusterNames

overwriteMetaclusterNames

## Usage

``` r
overwriteMetaclusterNames(fsom_parent, fsom_sub)
```

## Arguments

- fsom_parent:

  Table with all cells of interest

- fsom_sub:

  Subset of table used for reclustering

  The rows in `fsom_parent` and `fsom_sub` are matched to each other
  based on the `cell_id` column. Note that the numbers in the `cell_id`
  column of `fsom_sub` correspond to row indices of the root clustering.
  i.e. if `cell_id` is equal to 27 in `fsom_sub`, then that cell/row
  corresponds to row 27 in the clustering containing all subsets of
  cells (which unless you are doing multiple reclusterings, will most
  often be `fsom_parent`).

  Thus if `fsom_parent` does not have the column `cell_id`, a temporary
  one is made using its row indices.

## Value

A table with updated metacluster names from the reclustering.
