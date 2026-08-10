# gatingSetToTable

Utility function for converting between data types

## Usage

``` r
gatingSetToTable(gs, population = "root")
```

## Arguments

- gs:

  Preprocessed `GatingSet` to convert

- population:

  A string giving the name of the gated population to use for table
  creation

  Two columns are added to the returned `data.table` that are not
  present in the `GatingSet`:

  - `.id`, indicating which sample the cell belongs to

  - `cell_id`, indicating a unique key for each cell in the table.
    Mostly applicable when subsetting

## Value

A `data.table`, created from the given `GatingSet`
