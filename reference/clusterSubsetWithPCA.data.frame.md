# clusterSubsetWithPCA.data.frame

Recluster a subset of flow data using principal components of a data
frame.

## Usage

``` r
# S3 method for class 'data.frame'
clusterSubsetWithPCA(
  input,
  pca_obj,
  num_components,
  num_clus,
  xdim = 10,
  ydim = 10,
  seed = NULL,
  convert_to_channels = TRUE,
  fsom_file = NULL
)
```
