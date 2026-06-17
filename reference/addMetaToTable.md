# addMetaToTable

addMetaToTable

## Usage

``` r
addMetaToTable(table, sample_dt, join_col)
```

## Arguments

- table:

  Table containing expression data for all samples

- sample_dt:

  Table containing sample info

- join_col:

  Column on which to join the two tables described above. Or, to join on
  different variables between tables, an expression (e.g. `a == b`).

## Value

The data.table provided as input, with columns for sample information
added.
