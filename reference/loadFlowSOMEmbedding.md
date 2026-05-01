# loadFlowSOMEmbedding

Load saved FlowSOM and apply to new data

## Usage

``` r
loadFlowSOMEmbedding(file, newData)
```

## Arguments

- file:

  A file path to the FlowSOM object whose clusters you would like
  `newData` to map to

- newData:

  A flowFrame, flowSet, matrix with column names, or a list of paths for
  files or directories containing the data to be clustered

## Value

A FlowSOM object
