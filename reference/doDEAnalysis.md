# doDEAnalysis

doDEAnalysis

## Usage

``` r
doDEAnalysis(
  input,
  cols_to_test,
  design,
  contrasts,
  subpopulations = NULL,
  inverse = FALSE
)
```

## Arguments

- input:

  A `GatingSet` or `data.table`

- cols_to_test:

  A `character` vector specifying which channels should be tested for
  differential expression. Each element should be a column in the data's
  expression matrix.

- design:

  A design matrix from
  [`makeDesignMatrix()`](https://00berst33.github.io/flowFun/reference/makeDesignMatrix.md)

- contrasts:

  A contrasts matrix from
  [`makeContrastsMatrix()`](https://00berst33.github.io/flowFun/reference/makeContrastsMatrix.md)

- subpopulations:

  Only valid if `input` is a `GatingSet`. A `character` vector
  specifying which populations in the `GatingSet` to include in tests.

- inverse:

  A `boolean`, whether or not data should be transformed back to a
  linear scale before calculating MFIs and performing tests. Only valid
  if `input` is a `GatingSet` with a transformation and its inverse
  attached.

## Value

An `MArrayLM` object resulting from
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)
