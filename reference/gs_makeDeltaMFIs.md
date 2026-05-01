# gs_makeDeltaMFIs

gs_makeDeltaMFIs

## Usage

``` r
gs_makeDeltaMFIs(gs, control_gs, subpopulations, cols, metadata_col = NULL)
```

## Arguments

- gs:

  A GatingSet to find delta MFIs for

- control_gs:

  A GatingSet containing control samples

- subpopulations:

  The populations to perform the calculation for. Passed to
  [`gs_makeMFIMatrix()`](https://00berst33.github.io/flowFun/reference/gs_makeMFIMatrix.md)

- cols:

  The markers to find delta MFIs for

- metadata_col:

  The column in `pData(gs)` to use to match control and primary samples.
  If `NULL`, samples in both `GatingSets` will be assumed to be in the
  correct order
