# addMetadataToGatingSet

addMetadataToGatingSet

## Usage

``` r
addMetadataToGatingSet(gs, sample_df)
```

## Arguments

- gs:

  A `GatingSet`

- sample_df:

  A `data.frame`, where each row corresponds to a sample and each column
  corresponds to metadata about the samples, such as experimental group,
  or the filename of a corresponding control sample. Must contain a
  column for filename, called `File.Name`.
