# gs_makeMFIMatrix

gs_makeMFIMatrix

## Usage

``` r
gs_makeMFIMatrix(gs, cols, subpopulations, inverse = FALSE)
```

## Arguments

- gs:

  A `GatingSet`

- cols:

  A `character` vector of the columns in the expression matrix to
  calculate MFIs for

- subpopulations:

  A `character` vector of subpopulations in `gs`

- inverse:

  `boolean`, whether or not the data should be transformed back to its
  raw scale. If `TRUE`, a transformation with an associated inverse must
  be attached to the `GatingSet`

## Value

A data.frame of MFIs
