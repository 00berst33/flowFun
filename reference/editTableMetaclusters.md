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
fsom_table <- editTableMetaclusters(fsom_table,
                                    new_labels = c(
                                    "9" = "1",
                                    "2" = "1",
                                    "8" = "4",
                                    "23" = "15",
                                    "21" = "15",
                                    "20" = "15",
                                    "22" = "15",
                                    "7" = "6"),
                                    cluster_assignments = c("4" = "Undefined")
                                    )
#> Error: object 'fsom_table' not found
```
