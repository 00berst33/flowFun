# makeDesignMatrix

Generate a design matrix. (edit example)

## Usage

``` r
makeDesignMatrix(sample_df)
```

## Arguments

- sample_df:

  A data frame, generated either manually or by
  [`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md).

## Value

A matrix where each column is a group of interest, and each row is a
sample.

## Examples

``` r
filepath <- system.file("extdata", "sample_information.csv", package = "flowFun")
samples <- prepareSampleInfo(filepath, "Sample.Name", "File.Name", comparisons)
#> Warning: file("") only supports open = "w+" and open = "w+b": using the former
#> Error in read.table(file = file, header = header, sep = sep, quote = quote,     dec = dec, fill = fill, comment.char = comment.char, ...): no lines available in input

design <- makeDesignMatrix(samples)
#> Error: object 'samples' not found
```
