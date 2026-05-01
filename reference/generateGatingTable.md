# generateGatingTable

generateGatingTable

## Usage

``` r
generateGatingTable(gs, collapse_data = TRUE, ld_stain = NULL)
```

## Arguments

- gs:

  A `GatingSet`

- collapse_data:

  A `boolean`, indicating if samples should be collapsed into one.
  Default is `TRUE`

- ld_stain:

  If the experiment contains a live/dead stain, a `character` indicating
  which channel

## Value

A data.table, to be used for creation of a `gatingTemplate`

## See also

[`openCyto::gatingTemplate()`](https://rdrr.io/pkg/openCyto/man/gatingTemplate-class.html)
