# plotLabeled2DScatter

Plot a 2D scatterplot colored by metacluster.

## Usage

``` r
plotLabeled2DScatter(
  input,
  channelpair,
  clusters = NULL,
  metaclusters = NULL,
  labels = TRUE,
  colors = NULL
)
```

## Arguments

- input:

  A FlowSOM object or data table.

- channelpair:

  A vector of two channel names.

- clusters:

  A vector specifying clusters of interest.

- metaclusters:

  A vector specifying metaclusters of interest.

- labels:

  Logical, should labels be plotted for each cluster center. Default is
  `TRUE`.

- colors:

  Optional, a vector of colors for each metacluster.

## Value

A plot drawn with ggplot2.

## Details

It is only necessary to define one of the arguments `clusters` and
`metaclusters`, but the user may use both if they wish to plot a
particular subset of their data.

This function draws a 2D scatterplot with the given channels on each
axis. Cells are colored by metacluster, and if `labels = TRUE`, a label
is plotted at each cluster center. The number contained in these labels
indicates the cluster's number, and the color of the label indicates
which metacluster it belongs to.

## Examples

``` r
# not run

# Generate plot
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV615-P-A"),
                     metaclusters = c(17, 19))
#> Error: object 'fsom_dt' not found

# Generate same plot with no labels
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV615-P-A"),
                     metaclusters = c(17, 19),
                     labels = FALSE)
#> Error: object 'fsom_dt' not found

# Specify both metaclusters and clusters
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV563-A"),
                     clusters = c(66, 79, 91),
                     metaclusters = 16)
#> Error: object 'fsom_dt' not found

# Specifying only clusters still colors by metacluster
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV563-A"),
                     clusters = c(66, 79, 91))
#> Error: object 'fsom_dt' not found
```
