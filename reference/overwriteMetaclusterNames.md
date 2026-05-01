# overwriteMetaclusterNames

overwriteMetaclusterNames

## Usage

``` r
overwriteMetaclusterNames(fsom_dt, fsom_sub)
```

## Arguments

- fsom_dt:

  Table with all cells of interest

- fsom_sub:

  Subset of table used for reclustering

  Both tables should have a `cell_id` and `Metacluster` column.

## Value

A table with updated metacluster names from the reclustering.
