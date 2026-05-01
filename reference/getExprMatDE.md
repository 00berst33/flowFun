# getExprMatDE

Find expression matrix, where rows are a metacluster/marker pair, and
columns are samples. Helper for
[`doDEAnalysis()`](https://00berst33.github.io/flowFun/reference/doDEAnalysis.md)

## Usage

``` r
getExprMatDE(fsom_dt, cols_to_use, sample_col = .id)
```

## Arguments

- fsom_dt:

  A data.table with columns for markers/channels, a column `File`
  denoting the sample a cell is from, and a column `Metacluster`
  denoting the metacluster a cell belongs to

- cols_to_use:

  A character vector of channels of interest; these should be column
  names of `fsom_dt`

- sample_col:

  The column in `fsom_dt` corresponding to sample ID. Assumed to be
  `.id` by default

  The resulting table is passed onto
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) in this
  package's typical workflow, but it may also serve as a helpful summary
  of the data.

## Value

A data.table where columns are samples and rows are metacluster/marker
groups
