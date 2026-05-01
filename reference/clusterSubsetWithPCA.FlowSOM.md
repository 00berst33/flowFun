# clusterSubsetWithPCA.FlowSOM

Recluster a subset of flow data using principal components of a filtered
aggregate .fcs file.

## Usage

``` r
# S3 method for class 'FlowSOM'
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
