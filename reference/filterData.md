# filterData

filterData

## Usage

``` r
filterData(input, clusters = NULL, metaclusters = NULL)
```

## Arguments

- input:

  A data table or FlowSOM object.

- clusters:

  A character vector specifying which clusters to keep. Default is
  `NULL`.

- metaclusters:

  A character vector specifying which metaclusters to keep. Default is
  `NULL`.

## Value

A data.table or flowFrame, where only cells belonging to the given
clusters or metaclusters are included.
