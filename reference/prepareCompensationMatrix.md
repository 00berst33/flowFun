# prepareCompensationMatrix

Prepares a compensation matrix for use with the
[`compensate()`](https://rdrr.io/pkg/flowCore/man/compensation-class.html)
function.

## Usage

``` r
prepareCompensationMatrix(matrix, gs)
```

## Arguments

- matrix:

  A `matrix` or `data.frame`; the compensation matrix to prepare

- gs:

  The `GatingSet` whose column names the compensation matrix's columns
  should be matched to

## Details

For the `flowWorkspace` function
[`compensate()`](https://rdrr.io/pkg/flowCore/man/compensation-class.html)
to work, the column names (checked by calling
[`colnames()`](https://rdrr.io/r/base/colnames.html) on an object) of
the `GatingSet` and compensation matrix must match. Oftentimes, by
default, there will be minor name differences between these two objects
e.g. two corresponding columns may be named `BUV805-A` and
`BUV805-A :: CD4` in the `GatingSet` and compensation matrix,
respectively.

This function simplifies the process of matching and editing such column
names by using substrings. In the example above, the compensation matrix
column `BUV805-A :: CD4` would search for a column in the `GatingSet`
named either `CD4` or `BUV805-A`, and be renamed to the match it finds.

An additional check that columns are in the correct order is done by
ensuring that the compensation matrix's main diagonal contains only 1s.

## Examples

``` r
# reading in compensation matrix
  comp_mat <- system.file("extdata", "auto_compensation_matrix.csv", package = "flowFun")
  comp_mat <- read.csv(comp_mat,
                       check.names = FALSE) # prevents column names from being edited

 # edit column names and reorder columns
  # not run: assume gs is a GatingSet
  prepareCompensationMatrix(comp_mat, gs)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colnames': object 'gs' not found
```
