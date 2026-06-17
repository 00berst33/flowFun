# Function to match a character vector whose elements are either markers or channels to a parent vector whose elements contain both markers and the corresponding channel. i.e. match "CD3" to "CD3 \<Alexa Fluor 700-A\>".

Function to match a character vector whose elements are either markers
or channels to a parent vector whose elements contain both markers and
the corresponding channel. i.e. match "CD3" to "CD3 \<Alexa Fluor
700-A\>".

## Usage

``` r
getFullNames(substr, target_names, keep_names = FALSE)
```
