# plotUMAP

Plots UMAP colored by metacluster.

## Usage

``` r
plotUMAP(input, num_cells = 5000, labels = NULL, colors = NULL, seed = NULL)
```

## Arguments

- input:

  A FlowSOM object or data table.

- num_cells:

  The number of cells to use for the dimension reduction. Default is
  5000.

- labels:

  Optional, a vector specifying the order in which metacluster name will
  appear in the legend. All metacluster names present in `input` should
  exists in this vector.

- colors:

  Optional, a vector of colors to use in the plot for each metacluster.

- seed:

  Optional, a seed for reproducibility.

## Value

A UMAP plot drawn with ggplot2.
