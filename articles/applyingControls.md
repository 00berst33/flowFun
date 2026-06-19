# Applying controls to samples

If you have control samples, they may be used to calculate delta MFIs,
which can then be used to plotting and differential expression testing.
Doing so requires a few more steps than the basic ones outlined in the
main workflow article.

Before doing the steps we outline below, users should have a `GatingSet`
corresponding to the primary data which the controls will be applied to.
This `GatingSet` should be complete with gates created during
pre-processing and clustering via the methods outlined in the main
workflow article. Recall that `GatingSets` contain information including
transformations, compensations, and gating schemes; this feature will be
used to quickly process the data of the control files.

### Reading and preprocessing control files

Like in the main workflow, the control files are read in with the
function `load_cytoset_from_fcs()`. Instead of creating a fresh
`GatingSet` with the data, however, we use the function
`gh_apply_to_cs()`. This function takes the compensation matrices,
transformations, and gates decided when pre-processing the main data and
applies it to the control files, removing the need for the user to
remember and reapply these steps themselves.

``` r

# Specify path to control FCS files
ctrl_dir <- "C:/path/to/controls"
# Read in files as cytoset
ctrl <- flowWorkspace::load_cytoset_from_fcs(list.files(ctrl_dir, full.names = TRUE), which.lines = 20000)

# Load saved GatingSet corresponding to main data
gs1 <- flowWorkspace::load_gs("path/to/gatingset")
# Apply transformations, compensations and gates to control gs
ctrl_gs <- flowWorkspace::gh_apply_to_cs(gs1[[1]], ctrl, compensation_source = "template") # make sure to exclude boolean
```

The result is a pre-processed `GatingSet`, similar to the original one,
but now applied to the control files.

### Cluster controls

Before we can calculate delta MFIs, it is necessary to categorize the
cells of the control files into the cell types identified in the main
data. At this point there should be no need for the over-clustering and
merging done in the main workflow; we simply use the final clustering
from the main data and apply it to different data.

It is assumed that the main `GatingSet`, in this case `gs1`, contains a
reference to an RDS file with information about the clustering. This
should have been added with a call to `addClusteringToGatingSet()`; you
can check that this is the case with `pData(gs1)`.

Below we use the function
[`clusterControls()`](https://00berst33.github.io/flowFun/reference/clusterControls.md)
to cluster the control files, which returns a `FlowSOM` object. We
convert this to a `data.table` and save results with
[`addClustersToGatingSet()`](https://00berst33.github.io/flowFun/reference/addClustersToGatingSet.md).

``` r

# Apply clustering to controls
# This may also be done manually using `FlowSOM::NewData`
fsom_projected <- clusterControls(ctrl_gs, gs1, "live")

# Get data.table for controls, and add clusters to corresponding GatingSet
ctrl_dt <- flowSOMToTable(fsom_projected)
addClustersToGatingSet(ctrl_dt, ctrl_gs, "live")
```

### Adding metadata to samples

In order to apply the controls to the data, the pipeline needs
information about which sample each control file corresponds to. We
incorporate this information with the function
[`addMetadataToGatingSet()`](https://00berst33.github.io/flowFun/reference/addMetadataToGatingSet.md).
To use it, you must have a table with sample information (see the
vignette on preparing data for differential analysis for more details on
how to create this table).

The table must have at least a column `filename`, containing filenames
of each sample, and a column(s) specifying control files and which
sample file they correspond to. An example is shown below, where the
`"FMO"` column lists the control files.

| sample.name | filename     | treatment   | FMO        |
|:------------|:-------------|:------------|:-----------|
| Sample 1    | sample_1.fcs | Treatment 1 | T1_FMO.fcs |
| Sample 2    | sample_2.fcs | Treatment 2 | T2_FMO.fcs |
| Sample 3    | sample_3.fcs | Treatment 2 | T2_FMO.fcs |
| Sample 4    | sample_4.fcs | Treatment 1 | T1_FMO.fcs |
| Sample 5    | sample_5.fcs | Treatment 2 | T2_FMO.fcs |
| Sample 6    | sample_6.fcs | Treatment 1 | T1_FMO.fcs |

Below we call `addMetaDataToGatingSet()`, which uses this table to add
metadata to the `GatingSet`’s `pData()`.

Note: adding this information to `pData(gs1)` allows the use of
functions in the `ggcyto` package to facet by any group defined in the
metadata.

``` r

# Read in table of sample info
# Should have column named 'filename'
sample_info <- read.csv("path/to/sample/info/csv")
# Add metadata to GatingSet
addMetadataToGatingSet(gs1, sample_info)
# Check results (not run)
pData(gs1)
```

### Calculating delta MFIs

Finally, we can calculate delta MFIs with the function
[`gs_makeDeltaMFIs()`](https://00berst33.github.io/flowFun/reference/gs_makeDeltaMFIs.md).
It is necessary to specify which column of the metadata table
(`pData(gs1)`) contains the names of the control files to use for this
calculation. In this case it is done by setting `metadata_col = "FMO"`.

``` r

# Get delta MFIs
delta_mfis <- gs_makeDeltaMFIs(gs1,
                               ctrl_gs,
                               subpopulations = subpops,
                               metadata_col = "FMO") # the name of the column containing control filenames
                                                     # (must be present in pData(gs1))
```
