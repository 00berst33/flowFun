# flowFun

``` r

library(flowFun)
# this should be the "getting started" vignette
```

To improve the accessibility of the package, `flowFun` provides a number
of template scripts for users to make getting started less daunting.
When viewing the GitHub page, these may be found in the `\scripts`
folder.

`flowFun` was also designed to allow users to easily import and export
data at any stage of the analysis. This page will give a brief overview
of the pipeline and how users may import their data at each step. The
diagram below summarizes each function relevant to the package, and
which step they fit into.

The most important data type in this pipeline is the `GatingSet`, which
comes from the `flowWorkspace` package. Users who have performed flow
analysis in R before may be familiar with the `flowSet`, which is quite
similar.
