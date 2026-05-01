# makeContrastsMatrix

Generate a contrasts matrix. (edit example)

## Usage

``` r
makeContrastsMatrix(sample_df, comparisons)
```

## Arguments

- sample_df:

  A data frame, generated either manually or by
  [`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md).

- comparisons:

  A list of comparisons, where each comparison is further defined by a
  list of factor levels. Each comparison and group of factor levels
  should be named. See examples for further detail on how this parameter
  should be defined.

## Value

A matrix, where each column corresponds to a comparison, and each row
corresponds to a group.

## Examples

``` r
comparisons <- list(
  male_vs_female = list(Sex = list("male", "female")),
  male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female"))
)

filepath <- system.file("extdata", "sample_information.csv", package = "flowFun")
samples <- prepareSampleInfo(filepath, "Sample.Name", "File.Name", comparisons)
#> Warning: file("") only supports open = "w+" and open = "w+b": using the former
#> Error in read.table(file = file, header = header, sep = sep, quote = quote,     dec = dec, fill = fill, comment.char = comment.char, ...): no lines available in input

contrasts <- makeContrastsMatrix(samples, comparisons)
#> Error: object 'samples' not found
```
