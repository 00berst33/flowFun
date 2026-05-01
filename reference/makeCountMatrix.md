# makeCountMatrix

Generate matrix of sample/metacluster cell counts.

## Usage

``` r
makeCountMatrix(
  input,
  sample_df = NULL,
  meta_names = NULL,
  min_cells = 3,
  min_samples = NULL
)
```

## Arguments

- input:

  A FlowSOM object or data frame.

- sample_df:

  If `input` is a FlowSOM object, a data frame from
  [`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md).

- meta_names:

  A vector of metacluster names of interest. By default, all metacluster
  names are used.

- min_cells:

  An integer, the minimum number of cells a metacluster should have in a
  specified number of samples to be included in the analysis.

- min_samples:

  An integer, the minimum number of samples a metacluster should have at
  least `min_cells` events in to be included in the analysis. By
  default, this is half the total number of samples.

## Value

A matrix, where each column represents a sample, and each row represents
a metacluster.

## See also

[`makeCountMatrix.FlowSOM()`](https://00berst33.github.io/flowFun/reference/makeCountMatrix.FlowSOM.md)
