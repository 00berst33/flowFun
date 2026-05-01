# flagChannelNames

Checks for typos/discrepancies between a flowSet or cytoset's channel
names between samples

## Usage

``` r
flagChannelNames(input)
```

## Arguments

- input:

  A flowSet or cytoset

## Value

A list. For each sample, contains a matrix of channel/marker pairs that
were inconsistent, i.e. not found in all samples
