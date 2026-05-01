# plotPCAScree

Plot a scree plot for an aggregate .fcs file's principal components.

## Usage

``` r
plotPCAScree(pca_obj)
```

## Arguments

- pca_obj:

  A [`prcomp`](https://rdrr.io/r/stats/prcomp.html) object as generated
  by
  [`doPCA()`](https://00berst33.github.io/flowFun/reference/doPCA.md).

## Value

A scree plot generated with
[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html).

## Details

This function produces a scree plot, which may be used to determine
which principal components should be used in the analysis. Typically,
there is a distinct "elbow" in the plot, where the amount of variance
explained by each component becomes drastically smaller. So, if this
elbow appears at principal component 5, you would use components 1-5 to
recluster your data. If there is no obvious elbow, it is instead
reasonable to use whichever number of components explains roughly 80% of
the variance.
