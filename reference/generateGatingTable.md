# generateGatingTable

Generates a default gating template for preprocessing

## Usage

``` r
generateGatingTable(gs, collapse_data = TRUE, ld_stain = NULL)
```

## Arguments

- gs:

  A `GatingSet`

- collapse_data:

  A `boolean`, indicating if samples should be collapsed into one when
  drawing gates. Default is `TRUE`

- ld_stain:

  If the experiment contains a live/dead stain, a `character` indicating
  which channel

## Value

A `data.table`, to be used for creation of a `gatingTemplate`

## Details

The package `openCyto` provides functions to automatically draw gates
for cytometry data. To work, it requires a
[`gatingTemplate`](https://rdrr.io/pkg/openCyto/man/gatingTemplate-class.html)
object specifying basic information about the gates to be drawn.

This function generates a simple gating scheme for use with `openCyto`
containing the following steps:

- Remove margin events

- Remove debris

- Remove doublets

- Remove dead cells (if a viability stain was used in the panel)

The argument `collapse_data` affects how gates are drawn. If `TRUE`
(default), all samples are aggregated into one, and gates are drawn
using said aggregated sample. If `FALSE`, gates are instead drawn
separately for each sample.

## See also

[`gatingTemplate()`](https://rdrr.io/pkg/openCyto/man/gatingTemplate-class.html),
[`gt_gating()`](https://rdrr.io/pkg/openCyto/man/gatingTemplate-class.html)
