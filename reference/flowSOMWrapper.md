# flowSOMWrapper

flowSOMWrapper

## Usage

``` r
flowSOMWrapper(
  input,
  cols_to_cluster,
  num_clus,
  seed = NULL,
  fsom_file = NULL,
  ...
)
```

## Arguments

- input:

  A `flowSet`

- cols_to_cluster:

  A numeric or character vector specifying which columns should be used
  for clustering.

- num_clus:

  The number of metaclusters to create.

- seed:

  Optional, a seed for reproducibility. Default is `NULL`.

- fsom_file:

  An .rds filename. If not `NULL` (default), the resulting FlowSOM
  object will be saved under this name.

- ...:

  Additional parameters to pass to
  [`FlowSOM()`](https://rdrr.io/pkg/FlowSOM/man/FlowSOM.html).

## Value

A data.table annotated with cluster and metacluster assignments.

## Examples

``` r
# Define markers/columns to use for clustering
cols_to_cluster <- c("BB700-P-A", "APC-A", "Alexa Fluor 700-A", "APC-Cy7-A", "BV480-A",
                     "BV605-A", "BV650-A", "BV711-A", "BV750-P-A", "BV786-A", "BUV395-A",
                     "BUV563-A", "BUV615-P-A", "BUV661-A", "BUV737-A", "BUV805-A", "PE-A",
                     "PE-CF594-A", "PE-Cy5-A", "PE-Cy7-A")

# not run: assume `fs` is a previously created flowSet
fsom_dt <- flowSOMWrapper(fs,
                          cols_to_cluster = cols_to_cluster,
                          num_clus = 35, # create 35 metaclusters
                          xdim = 12,
                          ydim = 12, # use a 12x12 grid, i.e. 144 nodes in the SOM
                          seed = 42) # set a seed for reproducibility
#> Error: object 'fs' not found
```
