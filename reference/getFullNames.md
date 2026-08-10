# getFullNames

getFullNames

## Usage

``` r
getFullNames(substr, target_names, keep_names = FALSE)
```

## Arguments

- substr:

  A character vector of names to be matched

- target_names:

  A character vector of names that the elements of `substr` will attempt
  to be matched to. Should contain longer names than `substr`

- keep_names:

  Boolean, whether the returned vector should have names indicating
  which element of `substr` was matched to the full name

## Details

Function to match a character vector whose elements are either markers
or channels to a parent vector whose elements contain both markers and
the corresponding channel. i.e. match "CD3" to "CD3 \<Alexa Fluor
700-A\>".
