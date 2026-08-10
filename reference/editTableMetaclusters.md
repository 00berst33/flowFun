# editTableMetaclusters

Merge or relabel metaclusters, and reassign clusters.

## Usage

``` r
editTableMetaclusters(
  table,
  new_labels = NULL,
  cluster_assignments = NULL,
  level_order = NULL
)
```

## Arguments

- table:

  A data table with columns for "Metacluster" and "Cluster".

- new_labels:

  Optional. A named vector, where names are the original metacluster
  names, and values are the new labels.

- cluster_assignments:

  Optional. A named vector, where names are the cluster numbers, and
  values are metacluster names.

- level_order:

  A character vector giving the desired order of the metacluster labels.

## Value

A data.table with an added column `Meta_edited`.

## Examples

``` r
# not run: assume `fsom_table` is the result of a call to `flowSOMWrapper()`
fsom_table <- editTableMetaclusters(fsom_table,
                                    new_labels = c("9" = "1",
                                                   "2" = "1", # merge metaclusters 2 and 9 into 1
                                                   "8" = "4", # merge metacluster 8 and 4
                                                   "23" = "15",
                                                   "21" = "15",
                                                   "20" = "15",
                                                   "22" = "15", # merge metaclusters 23,21,20, and 22 into 15
                                                   "7" = "6"), # merge metaclusters 7 and 6
                                    cluster_assignments = c("87" = "Undefined")
                                    )
#> Error: object 'fsom_table' not found
```
