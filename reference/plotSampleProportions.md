# plotSampleProportions

Plot cell type proportions by sample.

## Usage

``` r
plotSampleProportions(count_mat)
```

## Arguments

- count_mat:

  A data frame or table, where rows are clusters and columns are
  samples.

## Value

A stacked bar plot made with
[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html).

## Details

By default, the names used for each sample will be its file name.

## See also

[`getExprMatDE()`](https://00berst33.github.io/flowFun/reference/getExprMatDE.md)
