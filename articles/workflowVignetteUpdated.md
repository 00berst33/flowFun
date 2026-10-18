# Main workflow

## Installation

Data used for this vignette is kept in a separate R package,
`flowFunData`. To install `flowFun` with built vignettes, users should
run the following code.

``` r

install.packages("devtools")

# Install with vignettes
#   The package 'flowFunData' must first be installed
devtools::install_github("00berst33/flowFunData")
#   Then install flowFun
devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = TRUE)
```

Otherwise, to get started with the package more swiftly, the user may
install `flowFun` without its vignettes.

``` r

install.packages("devtools")

# Install without vignettes
#   No installation of 'flowFunData', set argument `build_vignettes` to FALSE
devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = FALSE)
```

## Data Import

Flow Cytometry data often contains millions of cells, so to facilitate
fast and efficient manipulation of these large datasets, this package
implements data structures from the `flowWorkspace` package;
particularly the `cytoframe`, `cytoset`, `GatingHierarchy`, and
`GatingSet`.

The `cytoframe` and `cytoset` are analogous to the more well known
`flowSet` and `flowFrame` structures, but instead of storing data in
memory they store only a pointer to the underlying data (i.e. they are
reference classes). The size of this pointer is negligible compared to
the size of cytometry datasets, making operations with `cytosets` and
`cytoframes` far more memory efficient.

The `GatingHierarchy` and `GatingSet` build upon the `cytoset` and
`cytoframe`. They store not just the single cell data, but also samples,
groups, transformations, compensation matrices, and gates, conveniently
all in one object. Other advantages are that they facilitate automatic
gating in R, and the importing of workspaces from external software like
FlowJo and Cytobank. These capabilities will be discussed further below.

Understanding these data structures in detail is not absolutely required
to use this pipeline, but for those that would like to learn more, see
the relevant documentation in `flowWorkspace`. An introductory vignette
explaining the basics of these objects may be brought up in RStudio by
entering
[`vignette(package="flowWorkspace", "flowWorkspace-Introduction")`](https://bioconductor.org/packages/release/bioc/vignettes/flowWorkspace/inst/doc/flowWorkspace-Introduction.html)
into the console.

Below, we demonstrate how data may be imported into R from two possible
sources: a set of FCS files, or a FlowJo/Diva/Cytobank workspace.

#### Reading in FCS files

Most often, a user has a directory or set of directories containing FCS
files that they would like to import into R. The chunk below
demonstrates how to prepare a set of FCS files as a `GatingSet`.

``` r

library(flowWorkspace)
library(flowFun)
library(data.table)

# Specify path to .fcs file
dir <- system.file("extdata", "samples", package = "flowFunData")

# Load cytoset from .fcs files
cs <- flowWorkspace::load_cytoset_from_fcs(path=dir, pattern=".fcs", truncate_max_range = TRUE)

# Make GatingSet
gs <- flowWorkspace::GatingSet(cs)
```

First the user specifies the name of the folder containing their FCS
files, which are then loaded into a `cytoset` with the function
[`load_cytoset_from_fcs()`](https://rdrr.io/pkg/flowWorkspace/man/load_cytoset_from_fcs.html).
This `cytoset` is then used to construct a corresponding `GatingSet`;
note that both objects point to the same underlying data.

Due to their nature as reference classes, working with data in a
`GatingSet` and/or `cytoset` requires slightly more consideration than
typical R objects. Our variables `gs` and `cs` point to the same
underlying data, so any changes made to one will also be applied to the
other. This concept applies to any copy or subset of them.

To demonstrate this, we first examine the channel names of the samples
in `cs`.

``` r

# Get cytoset channel names
colnames(cs)
#>  [1] "FSC-A"             "FSC-H"             "FSC-W"            
#>  [4] "SSC-A"             "SSC-H"             "SSC-W"            
#>  [7] "FITC-A"            "BB630-A"           "BB660-P-A"        
#> [10] "BB700-P-A"         "BB790-P-A"         "APC-A"            
#> [13] "Alexa Fluor 700-A" "APC-Cy7-A"         "BV421-A"          
#> [16] "BV480-A"           "BV570-A"           "BV605-A"          
#> [19] "BV650-A"           "BV711-A"           "BV750-P-A"        
#> [22] "BV786-A"           "BUV395-A"          "BUV496-A"         
#> [25] "BUV563-A"          "BUV615-P-A"        "BUV661-A"         
#> [28] "BUV737-A"          "BUV805-A"          "PE-A"             
#> [31] "PE-CF594-A"        "PE-Cy5-A"          "PE-Cy5.5-A"       
#> [34] "PE-Cy7-A"          "Time"
```

Next we change the first channel name and observe what happens to the
column names of `gs`.

``` r

# Change name of first channel in cytoset
colnames(cs)[1] <- "EDIT"

# Print resulting GatingSet column names
colnames(gs)
#>  [1] "EDIT"              "FSC-H"             "FSC-W"            
#>  [4] "SSC-A"             "SSC-H"             "SSC-W"            
#>  [7] "FITC-A"            "BB630-A"           "BB660-P-A"        
#> [10] "BB700-P-A"         "BB790-P-A"         "APC-A"            
#> [13] "Alexa Fluor 700-A" "APC-Cy7-A"         "BV421-A"          
#> [16] "BV480-A"           "BV570-A"           "BV605-A"          
#> [19] "BV650-A"           "BV711-A"           "BV750-P-A"        
#> [22] "BV786-A"           "BUV395-A"          "BUV496-A"         
#> [25] "BUV563-A"          "BUV615-P-A"        "BUV661-A"         
#> [28] "BUV737-A"          "BUV805-A"          "PE-A"             
#> [31] "PE-CF594-A"        "PE-Cy5-A"          "PE-Cy5.5-A"       
#> [34] "PE-Cy7-A"          "Time"
```

We see that despite never directly writing code to change the channel
names in `gs`, the channel name `"FSC-A"` has also been changed to
`"EDIT"`.

Without being aware of this property of reference classes, a user could
unknowingly make edits to their underlying data while performing
analysis. However, this error is easily avoided by making a “deep copy”
with the function `gs_clone()`. The resulting `GatingSet` is identical
to the one it was made from, other than the critical difference that it
does not point to the same underlying data. Edits may be made to the
clone freely, without the potential of unintentionally overwriting
original data or unsaved changes.

Users should make a clone of their `GatingSet` as below before
proceeding with analysis.

``` r

# Clone GatingSet
gs1 <- flowWorkspace::gs_clone(gs)
```

#### Reading in FlowJo workspace

Users with data in FlowJo, Diva, or Cytobank may have interest in the R
package `CytoML`. Workspaces from these softwares may be imported into R
as a `GatingSet`, with transformations, compensations, and gating
schemes preserved. This gives users who may prefer to perform
preprocessing steps in one of these programs the ability to easily move
between R and their program of choice.

The code below demonstrates how a user may import a workspace from
FlowJo as a GatingSet. All that is required is the path to an XML file.

``` r

# Open FlowJo workspace in R from .xml file
flowjo_file <- "path/to/flowjo.xml"
ws <- CytoML::open_flowjo_xml(file)

# Make GatingSet from FlowJo workspace
gs <- CytoML::flowjo_to_gatingset(ws,
                                  path = "path/to/fcs_files")
```

A `GatingSet` may also be exported as a workspace with the use of
Docker; interested users should see the GitHub repository for `CytoML`.

## Preprocessing

Once setup files have been cleaned and a compensation matrix has been
generated, the next step is to preprocess the raw data. In order, the
steps this pipeline implements are as follows:

1.  Removing margin events (i.e. values that fall outside the range the
    cytometer should be able to pick up).
2.  Removing doublets.
3.  Removing debris.
4.  Compensating the data.
5.  Transforming the data.
6.  Removing dead cells.

Each these steps, other than acquiring a compensation matrix, may be
left entirely automated if desired. However, oftentimes parameter tuning
is necessary for optimal results. Furthermore, some samples’ data may be
distributed in a way that is markedly different from the average, making
automated gating ineffective. For this reason the user also has the
option of manually drawing gates, either by simply applying the sample
one to each sample, or retroactively selecting samples for which the
automated gating performed poorly and re-drawing gates as necessary.
Although it is implemented, it is recommended to avoid manual gating
where possible. Doing so improves the objectivity and reproducibility of
the analysis.

This pipeline primarily uses the package `openCyto` to implement
preprocessing. `openCyto` can create an automated gating pipeline using
a `gatingTemplate`, which is specified by a user-created spread sheet
with 10 columns: `alias`, `pop`, `parent`, `dims`, `gating_method`,
`gating_args`, `collapseDataForGating`, `groupBy`,
`preprocessing_method`, and `preprocessing_args`. To see the details on
how the values for these columns should be entered, see the
`openCytoVignette` in the `openCyto` package. `flowFun` generates a
default `gatingTemplate` that will draw the basic gates, so most users
only need to understand how to edit this template, and not concern
themselves with making one from scratch.

The function
[`generateGatingTable()`](https://00berst33.github.io/flowFun/reference/generateGatingTable.md)
creates this default. It returns a `data.table` so that the user may
make edits more easily before finalizing the `gatingTemplate`.

The code below generates a table with the argument `collapse_data` set
to `FALSE`, so once the gating scheme is applied, gates will be uniquely
determined for each sample.

``` r

library(openCyto)

# Specify L/D stain
ld_stain <- "BUV496-A"

# Get number of samples in experiment
num_samples <- length(gs1)

# Generate default template
gt_table <- generateGatingTable(gs1, collapse_data = TRUE, ld_stain)
gt_table
#>         alias    pop     parent        dims   gating_method
#>        <char> <char>     <char>      <char>          <char>
#> 1: nonMargins      +       root FSC-A,SSC-A        boundary
#> 2:  nonDebris      + nonMargins       FSC-A gate_mindensity
#> 3:   singlets      +  nonDebris FSC-A,FSC-H     singletGate
#> 4:       live      -   singlets    BUV496-A gate_mindensity
#>                        gating_args collapseDataForGating groupBy
#>                             <char>                <lgcl>   <int>
#> 1: min=c(0,0),max=c(262143,262143)                  TRUE       6
#> 2:                            <NA>                  TRUE       6
#> 3:                            <NA>                  TRUE       6
#> 4:                            <NA>                  TRUE       6
#>    preprocessing_method preprocessing_args
#>                  <lgcl>             <lgcl>
#> 1:                   NA                 NA
#> 2:                   NA                 NA
#> 3:                   NA                 NA
#> 4:                   NA                 NA
```

It may be helpful to edit `gating_args`. Again, more details can be
found in the `openCyto` vignettes, but it is most relevant to know how
to adjust the range that the gate can be drawn in. How to do so depends
on the function being used to determine the gate, which is indicated in
the `gating_method` column.

For instance, note the gate defined by the first row of `gt_table`,
whose alias is `nonMargins`. In the `gating_method` column, we see that
the value is `"boundary"`. The parameters that define the range for this
type of gate are `min` and `max`. We choose
`gating_args="min=c(0,0),max=c(262143,262143)"`, and as a result, the
gate cannot be drawn outside of the interval \[0, 262143\] for FSC-A or
SSC-A.

Before applying these gates to the data, the user should supply any
compensation matrices and choose whatever transformation they would like
to apply. Additionally, if your panel includes a viability stain, it is
necessary to specify the corresponding channel, so that dead cells may
be gated out accurately.

Below, the compensation matrix is read in, the transformation applied,
and the automated gating applied and visualized.

``` r

library(ggcyto)

# Make gatingTemplate from gt_table
gt <- openCyto::gatingTemplate(gt_table)

# Compensation matrix csv
comp_mat <- system.file("extdata", "compensation_matrix.csv", package = "flowFunData")
comp_mat <- read.csv(comp_mat, check.names = FALSE)

# Compensate
# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.
#   Check that that is the case here:
colnames(comp_mat)
#>  [1] "FITC-A"            "BB630-A"           "BB700-P-A"        
#>  [4] "APC-A"             "Alexa Fluor 700-A" "APC-Cy7-A"        
#>  [7] "BV421-A"           "BV480-A"           "BV605-A"          
#> [10] "BV650-A"           "BV711-A"           "BV750-P-A"        
#> [13] "BV786-A"           "BUV395-A"          "BUV496-A"         
#> [16] "BUV563-A"          "BUV615-P-A"        "BUV661-A"         
#> [19] "BUV737-A"          "BUV805-A"          "PE-A"             
#> [22] "PE-CF594-A"        "PE-Cy5-A"          "PE-Cy7-A"
colnames(gs1)
#>  [1] "FSC-A"             "FSC-H"             "FSC-W"            
#>  [4] "SSC-A"             "SSC-H"             "SSC-W"            
#>  [7] "FITC-A"            "BB630-A"           "BB660-P-A"        
#> [10] "BB700-P-A"         "BB790-P-A"         "APC-A"            
#> [13] "Alexa Fluor 700-A" "APC-Cy7-A"         "BV421-A"          
#> [16] "BV480-A"           "BV570-A"           "BV605-A"          
#> [19] "BV650-A"           "BV711-A"           "BV750-P-A"        
#> [22] "BV786-A"           "BUV395-A"          "BUV496-A"         
#> [25] "BUV563-A"          "BUV615-P-A"        "BUV661-A"         
#> [28] "BUV737-A"          "BUV805-A"          "PE-A"             
#> [31] "PE-CF594-A"        "PE-Cy5-A"          "PE-Cy5.5-A"       
#> [34] "PE-Cy7-A"          "Time"

compensate(gs1, comp_mat)
#> A GatingSet with 6 samples


# Transform
trans <- flowCore::estimateLogicle(gs1[[1]], channels = colnames(comp_mat))
flowWorkspace::transform(gs1, trans)
#> A GatingSet with 6 samples


# Apply gating scheme
gt_gating(gt, gs1)


## Check one gate for each sample
# nonDebris gate
plotAllSamples(gs1, "FSC-A", "SSC-A", "nonMargins", "nonDebris")
```

![](workflowVignetteUpdated_files/figure-html/preprocessing-1.png)

``` r


# singlets gate
plotAllSamples(gs1, "FSC-A", "FSC-H", "nonDebris", "singlets")
```

![](workflowVignetteUpdated_files/figure-html/preprocessing-2.png)

``` r


# live cell gate
plotAllSamples(gs1, !!enquo(ld_stain), "FSC-A", "singlets", "live")
```

![](workflowVignetteUpdated_files/figure-html/preprocessing-3.png)

Here we see that the debris gate for the third sample is poorly placed,
and multiple gates for live/dead cells need adjustment. The algorithm
used to place gates for these populations, \`gate_mindensity’, works by
fitting a kernel density estimator to the data and finding the two
largest peaks. When only one peak is found, or a sample has a more
unique distribution, it may not perform as well without some additional
arguments. These gates may be fixed with a few different approaches.

The first is to adjust the `gatingTemplate`. In many cases, restricting
the range in which the gate can be drawn resolves this kind of issue. To
set a minimum and maximum value for a gate drawn by `gate_mindensity`,
edit the argument `gate_range`. Note that making changes to the
`gatingTemplate` does not immediately affect the `GatingSet` it was
applied to. First, the old gate(s) and any associated children must be
deleted within the `GatingSet` with
[`flowWorkspace::gs_pop_remove()`](https://rdrr.io/pkg/flowWorkspace/man/gs_pop_add.html)
(this does not delete the gates defined in the `gatingTemplate`). The
new `gatingTemplate` can then be applied with
[`openCyto::gatingTemplate()`](https://rdrr.io/pkg/openCyto/man/gatingTemplate-class.html)
and the deleted gates will be redrawn according to the adjustments. The
code below demonstrates this with the `data.table` created earlier, and
plots the new gates generated after applying the new gating scheme.

``` r

# Edit arguments of gate_mindensity
gt_table[2, "gating_args"] <- "gate_range=c(0,80000)" # Set range of debris gate to [0,80000]

# Remove gate(s) that will be redrawn
# Here we delete "nonDebris", which also deletes its children gates (singlets, live) 
flowWorkspace::gs_pop_remove(gs1, "nonDebris")

# Generate new gating template
gt <- openCyto::gatingTemplate(gt_table)

# Apply to data
openCyto::gt_gating(gt, gs1)

# Plot new nonDebris gates
plotAllSamples(gs1, "FSC-A", "SSC-A", "nonMargins", "nonDebris")
```

![](workflowVignetteUpdated_files/figure-html/edit-gatingtable-1.png)

Another option is to draw gates with collapsed data, instead of
calculating them per sample. This approach is particularly helpful when
the data does not have many cells. It is also more straightforward,
since the same gate is applied to each sample. The code below
demonstrates this, and generates plots of the new gate for live cells.

``` r

# Get new gating table
gt_table <- generateGatingTable(gs1, collapse_data = TRUE, ld_stain)

# Remove gate(s) that will be redrawn
flowWorkspace::gs_pop_remove(gs1, "nonDebris")

# Generate new gating template
gt <- openCyto::gatingTemplate(gt_table)

# Apply to data
openCyto::gt_gating(gt, gs1)

# Plot new L/D gates
plotAllSamples(gs1, !!enquo(ld_stain), "FSC-A", "singlets", "live")
```

![](workflowVignetteUpdated_files/figure-html/collapsed-gt-1.png)

Finally, we may also adjust a gate by manually redrawing it. This can be
done with the function `editGateManual()`, which implements
`CytoExploreR::cyto_gate_draw()` for this pipeline. A window where the
user may draw a new gate on the sample(s) of interest is brought up, and
once the gate has been drawn, the given `GatingSet` is edited
accordingly. An example of how this function may be used is shown below,
where we choose to manually draw and replace the `"live_cells"` gate for
all samples with a rectangle gate. Other supported gate types are
polygon, ellipse, threshold, boundary, interval, quadrant, and web; see
`?CytoExploreR::cyto_gate_draw()` for more details.

``` r

### Edit gate for one sample:
# Check gate before:
autoplot(gs1[[3]])

# Example of a call to manually redraw gates; not run
editGateManual(gs1,
               node = "nonDebris", # name of population to redraw gates for
               gate_type = "rectangle", # type of gate to draw
               sample_ids = 1) # vector of sample indices to redraw gates for, default is NULL (all samples)
```

Once preprocessing results appear satisfactory, the user can either
export results to FlowJo using the `CytoML` package (note: this requires
Docker to be installed) or continue with this pipeline to cluster data.

It is recommended that the user saves their `GatingSet` before
continuing. Recall that the earlier preprocessing steps were performed
on a deep copy of the data to make edits more easily reversible, and
avoid accidentally changing the underlying data. At this point the
preprocessing steps are final, and so they should be saved so that
progress is not lost. To save a `GatingSet`, create a new folder in an
appropriate directory, then call
[`flowWorkspace::save_gs()`](https://rdrr.io/pkg/flowWorkspace/man/save_gs.html).
Note that more than one `GatingSet` may not be saved to the same
directory.

``` r

# Specify name of directory to save to
path <- file.path(getwd(), "gating_set")

# Make directory if it doesn't exist
dir.create(path)

# Save GatingSet
flowWorkspace::save_gs(gs1, path = path) # Set `path` to the name of your chosen directory
```

After saving the user may safely close their R session and return later,
if desired. To load a previously saved `GatingSet`, use
[`flowWorkspace::load_gs()`](https://rdrr.io/pkg/flowWorkspace/man/save_gs.html).
The only input necessary is the name of the folder the `GatingSet` was
saved to.

``` r

# Load saved GatingSet
flowWorkspace::load_gs(path)
```

## Clustering

This package is compatible with various clustering algorithms, but was
designed specifically with FlowSOM in mind due to its fast runtime and
tools for data visualization. Therefore, much focus will be given to
FlowSOM in the following section.

The strategy for identifying cell type populations outlined by this
workflow is not entirely automated, and instead uses a human-in-the-loop
approach. Such an approach aims to reduce the burden of analyzing such a
large volume of data, while still incorporating expert knowledge to
ensure that results are biologically meaningful. To accomplish this, the
data is initially overclustered. The expert then examines the result of
the clustering through various relevant plots, taking into account not
only marker expression but also cluster sizes and mathematical distance,
and uses their best judgement to merge clusters until an appropriate
final clustering has been reached and labelled.

It is also necessary for the user to specify which markers are to be
used for clustering. There are a few considerations to account for when
selecting these, the most important being that it is highly recommended
that any markers the user intends to test for differential expression
are NOT used for clustering. There is a number of reasons for this, but
in short, it generally makes little sense to test for differential
expression of a marker in a cluster when the cluster is defined by high
or low expression of said marker (i.e. testing for differential
expression of CD8 in CD8 T cells). Furthermore, parameters like side and
forward scatter should not be included, and if a viability stain was
used, it should also be excluded.

Below, the preprocessed data from earlier is clustered using the FlowSOM
algorithm via the function
[`flowSOMWrapper()`](https://00berst33.github.io/flowFun/reference/flowSOMWrapper.md).
This outputs the same table that is given to it, but with the columns
`Metacluster` and `Cluster` joined to the right, giving the results of
the clustering. Specifying a seed is optional, but strongly recommended
for the sake of reproducibility.

`FlowSOM()` does not take `GatingSet` data structures as input, so to
use this function users should subset their data to the appropriate
population with
[`gs_pop_get_data()`](https://rdrr.io/pkg/flowWorkspace/man/gh_pop_get_data.html)
and pass the result to
[`cytoset_to_flowSet()`](https://rdrr.io/pkg/flowWorkspace/man/convert.html).
This table may then be passed to
[`flowSOMWrapper()`](https://00berst33.github.io/flowFun/reference/flowSOMWrapper.md).

``` r

# Read in table
fs <- flowWorkspace::gs_pop_get_data(gs1, y = "live")
fs <- flowWorkspace::cytoset_to_flowSet(fs)

# Define markers/columns to use for clustering
cols_to_cluster <- c("BB700-P-A", "APC-A", "Alexa Fluor 700-A", "APC-Cy7-A", "BV480-A", 
                     "BV605-A", "BV650-A", "BV711-A", "BV750-P-A", "BV786-A", "BUV395-A", 
                     "BUV563-A", "BUV615-P-A", "BUV661-A", "BUV737-A", "BUV805-A", "PE-A",
                     "PE-CF594-A", "PE-Cy5-A", "PE-Cy7-A")

# Perform clustering
fsom_dt <- flowSOMWrapper(fs,
                          cols_to_cluster = cols_to_cluster,
                          num_clus = 35,
                          xdim = 12,
                          ydim = 12,
                          seed = 42)
```

| Metacluster | Cluster |
|:-----------:|:-------:|
|     12      |   66    |
|     34      |   138   |
|      5      |   23    |
|      8      |   35    |
|     18      |   85    |

A number of different plots may be generated by this package to aid the
user in merging metaclusters, but arguably the most important are
heatmaps and dimension reduction plots. Below, a heatmap of MFIs by
metacluster is generated.

``` r

# Generate heatmap
plotMetaclusterMFIs(fsom_dt)
```

![](workflowVignetteUpdated_files/figure-html/first-heatmap-1.png)

This heatmap is what the user should be using as their main reference
when merging metaclusters. The dendrogram on the left-hand side, which
displays a hierarchal clustering of the metaclusters, is of particular
interest, and this hierarchy may be followed to perform the clustering.
However, this dendrogram should not be followed blindly. First, the
metaclusters that are determined to be similar by the clustering may not
be of biological interest; some cell type markers may be of greater
significance than others, and this is something that the clustering does
not take into account. Second, dendrograms may be created with different
linkage methods (single-linkage, average-linkage, centroid-linkage,
etc.), and these different methods may result in slightly different
hierarchies. This is not to say that the dendrogram is useless, only
that it has its limitations.

The user should also consider cluster size when merging. If a
metacluster is rather small, and differs from its neighboring clusters
in irrelevant markers, then it is reasonable to merge it. However, if
the metacluster is small but distinct in markers of interest, it may be
left unmerged.

Finally, if a metacluster does not represent any particular cell type of
interest, it may be dropped from the analysis entirely. For example,
depending on how strict you were with your gating, there may be some
debris remaining in your preprocessed data. The clustering algorithm
will likely cluster these together, and you may decide to exclude them
from further analysis.

To further aid the decision making process, a UMAP colored by
metacluster is generated below.

``` r

# Generate UMAP
plotUMAP(fsom_dt, num_cells = 5000, seed = 42)
```

![](images/workflow_final/umap1.png)

This plot is most useful to check metaclusters’ relative sizes, and how
similar they are to one another. The closer clusters (and cells) are to
each other, the more similar they are to each other. However, the
opposite is not necessarily true.

By examining these plots, it starts to become clear which metaclusters
are more similar to one another, and what their phenotype may be. For
example, we see that metaclusters 1, 2, and 7 are all somewhat closely
linked by the dendrogram in the heatmap, and are rather close to one
another in the UMAP. They all have a high expression of the markers CD3
and CD8, but differ noticeably in markers PD1, CCR7, and CD45RA. It
seems that these metaclusters are each a subset of CD8 T cells, and it
is reasonable to merge them. Depending on the question being
investigated, the user may wish to instead keep them separate, and
perhaps specify exhausted and non-exhausted CD8 T cells, which could
also be reasonable.

It should be noted that if a user does desire to look at a particular
cell type in greater resolution, this workflow allows for backgating on
and reclustering cells belonging to a particular metacluster; therefore
identifying every sub-population of interest is not necessary at this
step. For now, this example will find an initial low-resolution
clustering.

The function
[`editTableMetaclusters()`](https://00berst33.github.io/flowFun/reference/editTableMetaclusters.md)
can be used to merge and rename metaclusters. A new column,
`Meta_original`, is added to the given table, specifying the
metaclusters that were found when the clustering algorithm was first
called. The existing `Metaclusters` column is edited to contain the new
assignments.

``` r

# Merge metaclusters
fsom_dt <- editTableMetaclusters(fsom_dt, 
                                 new_labels = c("7" = "1",
                                                "2" = "1",
                                                "5" = "3",
                                                "4" = "3",
                                                "5" = "3",
                                                "6" = "3",
                                                "13" = "10",
                                                "22" = "17"))
```

| Meta_original | Metacluster | Cluster |
|:-------------:|:-----------:|:-------:|
|      12       |     12      |   66    |
|      34       |     34      |   138   |
|       5       |      3      |   23    |
|       8       |      8      |   35    |
|      18       |     18      |   85    |

A new heatmap may be generated to reflect the new clustering, and help
decide whether any further merging should be performed.

``` r

# Generate new heatmap
plotMetaclusterMFIs(fsom_dt, cols_to_cluster)
```

![](workflowVignetteUpdated_files/figure-html/second-plots-1.png)

Examining this heatmap suggests that metacluster 9 may be a subset of
gdT cells, but the somewhat low expression of TCRgd makes it difficult
to determine this with much certainty. Sometimes, once the most obvious
merging has been done, UMAPs and MFI heatmaps alone are no longer
sufficient to make informed decisions. Two other notable functions are
provided to give further info about the clusters, the first being
[`plotClusterMFIs()`](https://00berst33.github.io/flowFun/reference/plotClusterMFIs.md).
It functions almost identically to
[`plotMetaclusterMFIs()`](https://00berst33.github.io/flowFun/reference/plotMetaclusterMFIs.md),
except that its rows are clusters, and it may be used to look at the
clusters within specific metacluster(s) of interest via the parameter
`metaclusters`. This is demonstrated below by calling the function on
metacluster 9.

``` r

# Generate cluster MFI heatmap
plotClusterMFIs(fsom_dt, cols_to_cluster, metaclusters = 9)
```

![](workflowVignetteUpdated_files/figure-html/cluster-heatmap-1.png)

It seems that clusters 55, 57, 58, 69, and 70 have a particularly low
expression of TCRgd compared to other nodes in this metacluster. This
may be further investigated with the second notable function,
[`plotLabeled2DScatter()`](https://00berst33.github.io/flowFun/reference/plotLabeled2DScatter.md).
It creates a 2D scatterplot colored by metacluster, where each cluster
center is labeled with its cluster’s number. Below, this function is
called on these five clusters.

``` r

# Generate 2D scatterplot
plotLabeled2DScatter(fsom_dt, 
                     channelpair = c("APC-Cy7-A", "BUV563-A"), 
                     clusters = c(55, 57, 58, 69, 70),
                     metaclusters = NULL)
#> [1] "CD3 <APC-Cy7-A>"  "TCRgd <BUV563-A>"
```

![](workflowVignetteUpdated_files/figure-html/2d-scatter-1.png)

Cluster 55 clearly does not have a high enough expression of TCRgd to
justify classifying it as gdT cells. Additionally, although it has
appropriately high expression of TCRgd, cluster 70 also has high
expression of multiple unrelated phenotyping markers (as seen in the
heatmap), suggesting that it may be a population of cells that were
stuck together. In the next merging step performed below, these clusters
are reassigned to a new metacluster named `"Undefined"`. It is
recommended to avoid reassigning clusters as much as possible, as doing
so increases the subjectivity of the analysis, but it is ultimately up
to the user to determine when this step is appropriate.

``` r

# Merge metaclusters
fsom_dt <- editTableMetaclusters(fsom_dt,
                                 new_labels = c("19" = "14",
                                                "20" = "14",
                                                "21" = "14",
                                                "18" = "14",
                                                "25" = "14",
                                                "16" = "8",
                                                "11" = "8",
                                                "23" = "15"),
                                 cluster_assignments = c("55" = "Undefined",
                                                         "70" = "Undefined"))
```

Again, new plots are generated.

``` r

# New heatmap
plotMetaclusterMFIs(fsom_dt, cols_to_cluster)

# New UMAP
plotUMAP(fsom_dt, num_cells = 5000, seed = 42)
```

![](images/workflow_final/umap2.png)

Finally, once the clustering is complete, the remaining metaclusters may
be labelled by setting each of them equal to the desired name, instead
of another metacluster to be merged into.

``` r

# Label metaclusters
fsom_dt <- editTableMetaclusters(fsom_dt, new_labels = c("14" = "Monocytes",
                                                         "24" = "cDC",
                                                         "10" = "B cells",
                                                         "12" = "NK T cells",
                                                         "17" = "NK cells",
                                                         "1" = "CD8 T cells",
                                                         "9" = "gdT cells",
                                                         "3" = "CD4 T cells",
                                                         "15" = "pDC",
                                                         "8" = "Undefined"))

# Generate final heatmap
plotMetaclusterMFIs(fsom_dt, cols_to_cluster)
```

![](workflowVignetteUpdated_files/figure-html/new-labels-1.png)

An annotated heatmap displaying the merged clusters may be created with
the function
[`annotateMFIHeatmap()`](https://00berst33.github.io/flowFun/reference/annotateMFIHeatmap.md).

``` r

# Generate annotated heatmap
annotateMFIHeatmap(fsom_dt, cols_to_cluster)
```

![](workflowVignetteUpdated_files/figure-html/annotate-heatmap-1.png)

#### Saving results

The clustering results may now be saved to the `GatingSet` created
during preprocessing, as a set of boolean gates. This step is
implemented by the function
[`addClustersToGatingSet()`](https://00berst33.github.io/flowFun/reference/addClustersToGatingSet.md),
which needs only three inputs: the `data.table` resulting from the
clustering step, the `GatingSet` whose expression data was used for
clustering, and the name of the population/gate in the `GatingSet` to
add these new gates to.

To see the gates that are currently in a `GatingSet`, the user can use
the
[`gs_get_pop_paths()`](https://rdrr.io/pkg/flowWorkspace/man/gs_get_pop_paths.html)
function.

``` r

# See all populations in GatingSet
flowWorkspace::gs_get_pop_paths(gs1)
#> [1] "root"                                "/nonMargins"                        
#> [3] "/nonMargins/nonDebris"               "/nonMargins/nonDebris/singlets"     
#> [5] "/nonMargins/nonDebris/singlets/live"
```

Note that the function `openCyto::plot()` used above to visualize the
`gatingTemplate` may also be used on a `GatingSet` or `GatingHierarchy`
object, and may be helpful in determining the appropriate parent
population.

Once the parent population has been determined, the new gates may be
added with
[`addClustersToGatingSet()`](https://00berst33.github.io/flowFun/reference/addClustersToGatingSet.md).
Because `GatingSet` is a reference class, only the underlying data is
edited and so there is no object returned from the function.

We visualize the results with `openCyto::plot()`:

``` r

# Add clusters to GatingSet
addClustersToGatingSet(fsom_dt, gs1, parent_gate = "live")

# Visualize new gating template
openCyto::plot(gs1)
```

![](workflowVignetteUpdated_files/figure-html/add-gates-to-gatingset-1.png)

Each cluster created earlier is now saved in the `GatingSet`. This is
convenient for several reasons. First, it allows for swift and easy
subsetting of the `GatingSet` by any combination of clusters. Any
function compatible with the `cytoset` or `flowSet` becomes available
for these subsets; in the below example the function
[`flowWorkspace::gs_pop_get_count_fast()`](https://rdrr.io/pkg/flowWorkspace/man/gs_pop_get_count_fast.html)
is used to get frequencies of CD4 T cells and CD8 T cells in each
sample.

``` r

# Get only CD8 and CD4 T cells in GatingSet
flowWorkspace::gs_pop_get_count_fast(gs1, 
                                     subpopulations = c("/nonMargins/nonDebris/singlets/live/CD8 T cells", 
                                                        "/nonMargins/nonDebris/singlets/live/CD4 T cells"), 
                                     statistic = "freq", 
                                     format = "long")
#>             name                                      Population
#>           <char>                                          <char>
#>  1: sample_1.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#>  2: sample_1.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#>  3: sample_2.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#>  4: sample_2.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#>  5: sample_3.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#>  6: sample_3.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#>  7: sample_4.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#>  8: sample_4.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#>  9: sample_5.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#> 10: sample_5.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#> 11: sample_6.fcs /nonMargins/nonDebris/singlets/live/CD8 T cells
#> 12: sample_6.fcs /nonMargins/nonDebris/singlets/live/CD4 T cells
#>                                  Parent Frequency ParentFrequency
#>                                  <char>     <num>           <num>
#>  1: /nonMargins/nonDebris/singlets/live    0.0138          0.3138
#>  2: /nonMargins/nonDebris/singlets/live    0.0386          0.3138
#>  3: /nonMargins/nonDebris/singlets/live    0.0106          0.4768
#>  4: /nonMargins/nonDebris/singlets/live    0.0422          0.4768
#>  5: /nonMargins/nonDebris/singlets/live    0.1074          0.7378
#>  6: /nonMargins/nonDebris/singlets/live    0.0730          0.7378
#>  7: /nonMargins/nonDebris/singlets/live    0.2244          0.5032
#>  8: /nonMargins/nonDebris/singlets/live    0.0670          0.5032
#>  9: /nonMargins/nonDebris/singlets/live    0.1610          0.5550
#> 10: /nonMargins/nonDebris/singlets/live    0.1154          0.5550
#> 11: /nonMargins/nonDebris/singlets/live    0.2104          0.5216
#> 12: /nonMargins/nonDebris/singlets/live    0.0672          0.5216
```

It is also simple to extract a subset of data using one of these gates,
and use it to write new FCS files. For example:

``` r

# Get cytoset for population of interest
cs_subset <- flowWorkspace::gs_pop_get_data(gs1, y = "NK cells")

# Make flowSet
fs <- flowWorkspace::cytoset_to_flowSet(cs_subset)

# Write new FCS files to disk
flowCore::write.flowSet(fs, outdir = getwd())
```

Recall that the `GatingSet` includes compensations, transformations,
sample information, metadata, and gating schemes in addition to the
expression data. The goal is to keep all relevant information about the
experiment and its analysis neatly contained within one lightweight
object, so that once it has been created, all that is necessary to
resume analysis is a call to
[`flowWorkspace::load_gs()`](https://rdrr.io/pkg/flowWorkspace/man/save_gs.html).

This feature is especially useful when performing a reclustering of the
data, or applying control stains, which will be demonstrated further
below.

#### Reclustering

This script also allows the user to perform backgating and reclustering
on any identified cell types of particular interest. This clustering may
be done using either a selected number of principal components obtained
from PCA, or a new subset of markers defined by the user. The former
approach is demonstrated below to further investigate the cells
identified as CD8 T cells.

``` r

# Get table of only CD8 T cells
new_table <- createFilteredAggregate(fsom_dt, 
                                     num_cells = Inf, 
                                     metaclusters = "CD8 T cells",
                                     clusters = NULL)
```

One an appropriate subset has been created, the next step is to create a
scree plot, which may be used to determine which principal components
should be used in the analysis. Typically, there is a distinct “elbow”
in the plot, where the amount of variance explained by each component
becomes drastically smaller. So, if this elbow appears at principal
component 5, components 1-5 would be selected to recluster the data. If
there is no obvious elbow, it is instead reasonable to use whichever
number of components explains roughly 70-80% of the variance. First,
[`doPCA()`](https://00berst33.github.io/flowFun/reference/doPCA.md) is
called to perform PCA. Then,
[`plotPCAScree()`](https://00berst33.github.io/flowFun/reference/plotPCAScree.md)
creates the scree plot.

``` r

# Perform PCA
pca_obj <- doPCA(new_table, cols_to_cluster)

# Draw scree plot
plotPCAScree(pca_obj)
```

![](workflowVignetteUpdated_files/figure-html/pca-1.png)

Here, there appears to be a fairly distinct elbow at `M = 6`. To
recluster using these principal components, the following call to
[`clusterSubsetWithPCA()`](https://00berst33.github.io/flowFun/reference/clusterSubsetWithPCA.md)
is made, with `num_components = 6`. Fewer metaclusters are created than
in the initial clustering, as fewer cell types are expected.

``` r

# Recluster CD8 T cells
fsom_sub <- clusterSubsetWithPCA(new_table, 
                                 pca_obj = pca_obj, 
                                 num_components = 6,
                                 num_clus = 15,
                                 seed = 33)
```

Note that the plot functions demonstrated earlier function the same for
any subsets created. Since the current cells of interest are CD8 T
cells, it is appropriate to define a new subset of markers to create the
heatmap.

``` r

# Define new columns to use
# cols_of_interest <- c(15, 16, 18, 23:25, 28, 32)
cols_of_interest <- c("APC-Cy7-A", "Alexa Fluor 700-A", "BUV395-A", "BV750-P-A", "PE-A", "BV786-A", "BV480-A", "BUV615-P-A")

# Generate heatmap for subsetted cells
plotMetaclusterMFIs(fsom_sub, cols_of_interest)

# Likewise, generate UMAP
plotUMAP(fsom_sub, num_cells = 5000, seed = 33)
```

![](images/workflow_final/reclustered_heatmap.png)![](images/workflow_final/reclustered_umap.png)

As with the parent population above, metaclusters are merged.

``` r

fsom_sub <- editTableMetaclusters(fsom_sub, new_labels = c("1" = "CD45RA+, CCR7-",
                                                           "2" = "CD45RA+, CCR7-",
                                                           "3" = "CD45RA+, CCR7-",
                                                           "4" = "CD45RA+, CCR7-",
                                                           "5" = "CD45RA+, CCR7-",
                                                           "6" = "CD45RA+, CCR7-",
                                                           "7" = "CD45RA-, CCR7-",
                                                           "8" = "CD45RA+, CCR7-",
                                                           "9" = "CD45RA+, CCR7-",
                                                           "10" = "CD45RA+, CCR7+",
                                                           "11" = "CD45RA-, CCR7-",
                                                           "12" = "CD45RA-, CCR7-",
                                                           "13" = "CD45RA-, CCR7+",
                                                           "14" = "CD45RA+, CCR7+",
                                                           "15" = "CD45RA+, CCR7+"))

# Generate heatmap for subsetted cells
plotMetaclusterMFIs(fsom_sub, cols_of_interest)
```

![](workflowVignetteUpdated_files/figure-html/reclustering-edit-1.png)

If desired, the clusters found in this subset may overwrite their parent
cluster in the `data.table` used for the initial clustering. Meaning, in
this case, any cells previously labeled as `"CD8 T cells"` will instead
be labelled with the more specific cell type found in the reclustering
above. New heatmaps and dimension reduction plots with the updated
labels may then be generated.

The function
[`overwriteMetaclusterNames()`](https://00berst33.github.io/flowFun/reference/overwriteMetaclusterNames.md)
performs this operation. The code chunk below demonstrates this and
generates a new heatmap.

``` r

# Overwrite labels for CD8 T cells with those found in reclustering
fsom_over <- overwriteMetaclusterNames(fsom_dt, fsom_sub)

# Generate new UMAP
# should minimum cell count per cluster be enforced?
plotUMAP(fsom_over, num_cells = 5000, seed = 42)
```

![](images/workflow_final/umap4.png)

The finalized populations may also be saved to the `GatingSet`. In this
case, we have a reclustering of the CD8 T cell population, so we call
[`addClustersToGatingSet()`](https://00berst33.github.io/flowFun/reference/addClustersToGatingSet.md)
with `parent_gate = "CD8 T cells"`.

``` r

# Add clusters to GatingSet
addClustersToGatingSet(fsom_sub, gs1, parent_gate = "CD8 T cells")

# Visualize new gating template
openCyto::plot(gs1)
```

![](workflowVignetteUpdated_files/figure-html/save-reclustering-to-gatingset-1.png)

The resulting plot shows each cluster defined so far, and any
parent-child relationships.

## Differential Analysis

#### Setup

The final analysis step in this workflow is differential analysis.
Differential marker expression is implemented using limma, and
differential count analysis using edgeR, both R packages found on
Bioconductor.

Some preparation is necessary before any tests can be performed. First,
it is necessary to specify information about each sample, including any
FMO or isotype controls, and which groups of interest it belongs to
(e.g. control vs. treatment). This information may be provided as a .csv
file or `data.frame` to the function
[`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md),
which will prepare the table for further use by the workflow and ensure
that its values are R-friendly (no special characters like @ or \#, no
empty values, etc.).

Further explanation may be found in the “Setting up differential
analysis” article, but in short it is most important that the table has
at least these two columns, with the same names:

- `sample.name`, which specifies the name of each sample
- `filename`, which specifies the filename of each sample. These names
  should appear as they do when calling
  [`sampleNames()`](https://rdrr.io/pkg/Biobase/man/featureNames.html)
  on your `GatingSet` or `flowSet`.

Below is an example of what an appropriate .csv file might look like.

|   filename   | sample.name | disease |  NAC   |       FMO        |     Isotype      |
|:------------:|:-----------:|:-------:|:------:|:----------------:|:----------------:|
| sample_1.fcs |   Ctrl 1    |  Ctrl   |        | FMO_sample_1.fcs | Iso_sample_1.fcs |
| sample_2.fcs |   Ctrl 2    |  Ctrl   |        | FMO_sample_2.fcs | Iso_sample_2.fcs |
| sample_3.fcs |   MIBC 3    |  MIBC   | No NAC | FMO_sample_3.fcs | Iso_sample_3.fcs |
| sample_4.fcs |   MIBC 4    |  MIBC   | No NAC | FMO_sample_4.fcs | Iso_sample_4.fcs |
| sample_5.fcs |   MIBC 5    |  MIBC   |  NAC   | FMO_sample_5.fcs | Iso_sample_5.fcs |
| sample_6.fcs |   MIBC 6    |  MIBC   |  NAC   | FMO_sample_6.fcs | Iso_sample_6.fcs |

Before calling
[`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md),
the user must also specify the comparisons they wish to make between
groups, as a series of nested lists. Each nested list should include at
least two levels of a factor to be used for the comparison. Each
comparison may be named however the user desires, given that it is
R-friendly. Below, two comparisons, `ctrl_vs_mibc` and `nac_vs_no_nac`
are defined.

``` r

comparisons <- list(
  ctrl_vs_mibc = list(disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(disease = "MIBC", NAC = list("NAC", "No.NAC"))
)
```

The factor names (which would be “disease” and “NAC” in the above
example), unlike the names of comparisons, should be identical to their
column names in the .csv file or data frame given to
[`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md).
Similarly, the factor levels (which would be “Ctrl”, “MIBC”, “NAC”, and
“No NAC” in the above example), should be listed identical to how they
appear in the data frame. Note how these values match those present in
the table just shown.

Once these two items have been specified, the remainder of analysis is
rather straightforward.
[`prepareSampleInfo()`](https://00berst33.github.io/flowFun/reference/prepareSampleInfo.md)
is called, formatting the given table as necessary, and adding a new
column called `group`.

``` r

# Get filepath
file <- system.file("extdata", "sample_info.csv", package = "flowFunData")

# Prepare metadata for further analysis
sample_info <- prepareSampleInfo(file, 
                                 name_col = "sample.name",
                                 filename_col = "filename",
                                 comparisons = comparisons)
```

| sample.name | filename | disease | NAC | FMO | Isotype | group |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Ctrl 1 | sample_1.fcs | Ctrl | X | FMO_sample_1.fcs | Iso_sample_1.fcs | Ctrl_X |
| Ctrl 2 | sample_2.fcs | Ctrl | X | FMO_sample_2.fcs | Iso_sample_2.fcs | Ctrl_X |
| MIBC 3 | sample_3.fcs | MIBC | No.NAC | FMO_sample_3.fcs | Iso_sample_3.fcs | MIBC_No.NAC |
| MIBC 4 | sample_4.fcs | MIBC | No.NAC | FMO_sample_4.fcs | Iso_sample_4.fcs | MIBC_No.NAC |
| MIBC 5 | sample_5.fcs | MIBC | NAC | FMO_sample_5.fcs | Iso_sample_5.fcs | MIBC_NAC |
| MIBC 6 | sample_6.fcs | MIBC | NAC | FMO_sample_6.fcs | Iso_sample_6.fcs | MIBC_NAC |

Next, the design, contrasts, and count matrices are created. The current
population of interest (in this case, `fsom_dt`) should be given to
[`makeCountMatrix()`](https://00berst33.github.io/flowFun/reference/makeCountMatrix.md).
Furthermore, the parameter `meta_names` should be used to specify the
metaclusters to use for analysis. By default, all metaclusters are
included, but if the final merging included an undefined metacluster, or
some are of no interest, it is useful to specify this parameter.

``` r

# Generate design matrix
design <- makeDesignMatrix(sample_info)

# Print result
data.frame(design)
#>   Ctrl_X MIBC_NAC MIBC_No.NAC
#> 1      1        0           0
#> 2      1        0           0
#> 3      0        0           1
#> 4      0        0           1
#> 5      0        1           0
#> 6      0        1           0
```

``` r

# Generate contrasts matrix
contrasts <- makeContrastsMatrix(sample_info, comparisons)
# Print
contrasts
#>              Contrasts
#> Levels        MIBC_No.NAC_OR_MIBC_NAC__vs__Ctrl_X MIBC_NAC__vs__MIBC_No.NAC
#>   Ctrl_X                                     -1.0                         0
#>   MIBC_NAC                                    0.5                         1
#>   MIBC_No.NAC                                 0.5                        -1
```

``` r

# Generate matrix of sample/metacluster cell counts
counts <- makeCountMatrix(fsom_dt, 
                          min_cells = 3, 
                          min_samples = 4)
```

|  | sample_1.fcs | sample_2.fcs | sample_3.fcs | sample_4.fcs | sample_5.fcs | sample_6.fcs |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|
| CD8 T cells | 69 | 53 | 537 | 1122 | 805 | 1052 |
| B cells | 36 | 26 | 430 | 49 | 73 | 61 |
| NK T cells | 127 | 35 | 74 | 40 | 23 | 20 |
| Monocytes | 368 | 812 | 912 | 291 | 786 | 377 |
| pDC | 123 | 43 | 304 | 6 | 18 | 14 |
| NK cells | 57 | 92 | 197 | 60 | 59 | 25 |
| cDC | 63 | 8 | 8 | 1 | 3 | 1 |
| 26 | 9 | 4 | 13 | 1 | 7 | 8 |
| 27 | 166 | 234 | 124 | 41 | 74 | 117 |
| 28 | 21 | 50 | 35 | 22 | 45 | 57 |
| 29 | 105 | 364 | 109 | 254 | 76 | 129 |
| CD4 T cells | 193 | 211 | 365 | 335 | 577 | 336 |
| 30 | 10 | 19 | 8 | 10 | 25 | 33 |
| 31 | 51 | 221 | 53 | 124 | 47 | 71 |
| 32 | 10 | 13 | 10 | 14 | 10 | 10 |
| 33 | 24 | 19 | 24 | 29 | 20 | 27 |
| 34 | 8 | 31 | 10 | 16 | 8 | 13 |
| 35 | 10 | 42 | 12 | 36 | 15 | 28 |
| Undefined | 115 | 106 | 372 | 58 | 100 | 222 |
| gdT cells | 4 | 1 | 92 | 7 | 4 | 7 |

To instead get counts from a `GatingSet`, the function
[`flowWorkspace::gs_pop_get_count_fast()`](https://rdrr.io/pkg/flowWorkspace/man/gs_pop_get_count_fast.html)
with argument `format` set to `"wide"` can also quickly obtain a count
matrix in the necessary format.

``` r

# Get counts from GatingSet
gs_counts <- flowWorkspace::gs_pop_get_count_fast(gs1, format = "wide") 
gs_counts
#>                                                                sample_1.fcs
#> /nonMargins                                                            4971
#> /nonMargins/nonDebris                                                  2338
#> /nonMargins/nonDebris/singlets                                         2138
#> /nonMargins/nonDebris/singlets/live                                    1569
#> /nonMargins/nonDebris/singlets/live/26                                    9
#> /nonMargins/nonDebris/singlets/live/27                                  166
#> /nonMargins/nonDebris/singlets/live/28                                   21
#> /nonMargins/nonDebris/singlets/live/29                                  105
#> /nonMargins/nonDebris/singlets/live/30                                   10
#> /nonMargins/nonDebris/singlets/live/31                                   51
#> /nonMargins/nonDebris/singlets/live/32                                   10
#> /nonMargins/nonDebris/singlets/live/33                                   24
#> /nonMargins/nonDebris/singlets/live/34                                    8
#> /nonMargins/nonDebris/singlets/live/35                                   10
#> /nonMargins/nonDebris/singlets/live/B cells                              36
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         193
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                          69
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+            7
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-           43
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+            4
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-           15
#> /nonMargins/nonDebris/singlets/live/Monocytes                           368
#> /nonMargins/nonDebris/singlets/live/NK T cells                          127
#> /nonMargins/nonDebris/singlets/live/NK cells                             57
#> /nonMargins/nonDebris/singlets/live/Undefined                           115
#> /nonMargins/nonDebris/singlets/live/cDC                                  63
#> /nonMargins/nonDebris/singlets/live/gdT cells                             4
#> /nonMargins/nonDebris/singlets/live/pDC                                 123
#> root                                                                   5000
#>                                                                sample_2.fcs
#> /nonMargins                                                            4971
#> /nonMargins/nonDebris                                                  3411
#> /nonMargins/nonDebris/singlets                                         3164
#> /nonMargins/nonDebris/singlets/live                                    2384
#> /nonMargins/nonDebris/singlets/live/26                                    4
#> /nonMargins/nonDebris/singlets/live/27                                  234
#> /nonMargins/nonDebris/singlets/live/28                                   50
#> /nonMargins/nonDebris/singlets/live/29                                  364
#> /nonMargins/nonDebris/singlets/live/30                                   19
#> /nonMargins/nonDebris/singlets/live/31                                  221
#> /nonMargins/nonDebris/singlets/live/32                                   13
#> /nonMargins/nonDebris/singlets/live/33                                   19
#> /nonMargins/nonDebris/singlets/live/34                                   31
#> /nonMargins/nonDebris/singlets/live/35                                   42
#> /nonMargins/nonDebris/singlets/live/B cells                              26
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         211
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                          53
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+            8
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-           29
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+            1
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-           15
#> /nonMargins/nonDebris/singlets/live/Monocytes                           812
#> /nonMargins/nonDebris/singlets/live/NK T cells                           35
#> /nonMargins/nonDebris/singlets/live/NK cells                             92
#> /nonMargins/nonDebris/singlets/live/Undefined                           106
#> /nonMargins/nonDebris/singlets/live/cDC                                   8
#> /nonMargins/nonDebris/singlets/live/gdT cells                             1
#> /nonMargins/nonDebris/singlets/live/pDC                                  43
#> root                                                                   5000
#>                                                                sample_3.fcs
#> /nonMargins                                                            4943
#> /nonMargins/nonDebris                                                  4438
#> /nonMargins/nonDebris/singlets                                         4230
#> /nonMargins/nonDebris/singlets/live                                    3689
#> /nonMargins/nonDebris/singlets/live/26                                   13
#> /nonMargins/nonDebris/singlets/live/27                                  124
#> /nonMargins/nonDebris/singlets/live/28                                   35
#> /nonMargins/nonDebris/singlets/live/29                                  109
#> /nonMargins/nonDebris/singlets/live/30                                    8
#> /nonMargins/nonDebris/singlets/live/31                                   53
#> /nonMargins/nonDebris/singlets/live/32                                   10
#> /nonMargins/nonDebris/singlets/live/33                                   24
#> /nonMargins/nonDebris/singlets/live/34                                   10
#> /nonMargins/nonDebris/singlets/live/35                                   12
#> /nonMargins/nonDebris/singlets/live/B cells                             430
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         365
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                         537
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+           98
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-          244
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+           40
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-          155
#> /nonMargins/nonDebris/singlets/live/Monocytes                           912
#> /nonMargins/nonDebris/singlets/live/NK T cells                           74
#> /nonMargins/nonDebris/singlets/live/NK cells                            197
#> /nonMargins/nonDebris/singlets/live/Undefined                           372
#> /nonMargins/nonDebris/singlets/live/cDC                                   8
#> /nonMargins/nonDebris/singlets/live/gdT cells                            92
#> /nonMargins/nonDebris/singlets/live/pDC                                 304
#> root                                                                   5000
#>                                                                sample_4.fcs
#> /nonMargins                                                            4961
#> /nonMargins/nonDebris                                                  3077
#> /nonMargins/nonDebris/singlets                                         2770
#> /nonMargins/nonDebris/singlets/live                                    2516
#> /nonMargins/nonDebris/singlets/live/26                                    1
#> /nonMargins/nonDebris/singlets/live/27                                   41
#> /nonMargins/nonDebris/singlets/live/28                                   22
#> /nonMargins/nonDebris/singlets/live/29                                  254
#> /nonMargins/nonDebris/singlets/live/30                                   10
#> /nonMargins/nonDebris/singlets/live/31                                  124
#> /nonMargins/nonDebris/singlets/live/32                                   14
#> /nonMargins/nonDebris/singlets/live/33                                   29
#> /nonMargins/nonDebris/singlets/live/34                                   16
#> /nonMargins/nonDebris/singlets/live/35                                   36
#> /nonMargins/nonDebris/singlets/live/B cells                              49
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         335
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                        1122
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+          124
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-          390
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+          239
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-          369
#> /nonMargins/nonDebris/singlets/live/Monocytes                           291
#> /nonMargins/nonDebris/singlets/live/NK T cells                           40
#> /nonMargins/nonDebris/singlets/live/NK cells                             60
#> /nonMargins/nonDebris/singlets/live/Undefined                            58
#> /nonMargins/nonDebris/singlets/live/cDC                                   1
#> /nonMargins/nonDebris/singlets/live/gdT cells                             7
#> /nonMargins/nonDebris/singlets/live/pDC                                   6
#> root                                                                   5000
#>                                                                sample_5.fcs
#> /nonMargins                                                            4896
#> /nonMargins/nonDebris                                                  3791
#> /nonMargins/nonDebris/singlets                                         3513
#> /nonMargins/nonDebris/singlets/live                                    2775
#> /nonMargins/nonDebris/singlets/live/26                                    7
#> /nonMargins/nonDebris/singlets/live/27                                   74
#> /nonMargins/nonDebris/singlets/live/28                                   45
#> /nonMargins/nonDebris/singlets/live/29                                   76
#> /nonMargins/nonDebris/singlets/live/30                                   25
#> /nonMargins/nonDebris/singlets/live/31                                   47
#> /nonMargins/nonDebris/singlets/live/32                                   10
#> /nonMargins/nonDebris/singlets/live/33                                   20
#> /nonMargins/nonDebris/singlets/live/34                                    8
#> /nonMargins/nonDebris/singlets/live/35                                   15
#> /nonMargins/nonDebris/singlets/live/B cells                              73
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         577
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                         805
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+           29
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-          389
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+           88
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-          299
#> /nonMargins/nonDebris/singlets/live/Monocytes                           786
#> /nonMargins/nonDebris/singlets/live/NK T cells                           23
#> /nonMargins/nonDebris/singlets/live/NK cells                             59
#> /nonMargins/nonDebris/singlets/live/Undefined                           100
#> /nonMargins/nonDebris/singlets/live/cDC                                   3
#> /nonMargins/nonDebris/singlets/live/gdT cells                             4
#> /nonMargins/nonDebris/singlets/live/pDC                                  18
#> root                                                                   5000
#>                                                                sample_6.fcs
#> /nonMargins                                                            4922
#> /nonMargins/nonDebris                                                  3648
#> /nonMargins/nonDebris/singlets                                         3464
#> /nonMargins/nonDebris/singlets/live                                    2608
#> /nonMargins/nonDebris/singlets/live/26                                    8
#> /nonMargins/nonDebris/singlets/live/27                                  117
#> /nonMargins/nonDebris/singlets/live/28                                   57
#> /nonMargins/nonDebris/singlets/live/29                                  129
#> /nonMargins/nonDebris/singlets/live/30                                   33
#> /nonMargins/nonDebris/singlets/live/31                                   71
#> /nonMargins/nonDebris/singlets/live/32                                   10
#> /nonMargins/nonDebris/singlets/live/33                                   27
#> /nonMargins/nonDebris/singlets/live/34                                   13
#> /nonMargins/nonDebris/singlets/live/35                                   28
#> /nonMargins/nonDebris/singlets/live/B cells                              61
#> /nonMargins/nonDebris/singlets/live/CD4 T cells                         336
#> /nonMargins/nonDebris/singlets/live/CD8 T cells                        1052
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7+           63
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA+, CCR7-          291
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7+          169
#> /nonMargins/nonDebris/singlets/live/CD8 T cells/CD45RA-, CCR7-          529
#> /nonMargins/nonDebris/singlets/live/Monocytes                           377
#> /nonMargins/nonDebris/singlets/live/NK T cells                           20
#> /nonMargins/nonDebris/singlets/live/NK cells                             25
#> /nonMargins/nonDebris/singlets/live/Undefined                           222
#> /nonMargins/nonDebris/singlets/live/cDC                                   1
#> /nonMargins/nonDebris/singlets/live/gdT cells                             7
#> /nonMargins/nonDebris/singlets/live/pDC                                  14
#> root                                                                   5000
```

By default, counts for all populations in the gating scheme are found.
However, for differential analysis, only the populations to be tested
should be used. Recall that the populations in a `GatingSet` may be
checked with
[`flowWorkspace::gs_get_pop_paths()`](https://rdrr.io/pkg/flowWorkspace/man/gs_get_pop_paths.html),
and that the full gating scheme may be visualized with
`openCyto::plot()`. The function
[`flowWorkspace::gs_pop_get_children()`](https://rdrr.io/pkg/flowWorkspace/man/gs_pop_get_children.html)
is used below to retrieve the paths for each cell type found in the
first round of clustering, as they all share the same parent. The count
matrix is then subsetted to the appropriate populations.

``` r

# Get paths to cluster populations
subpops <- flowWorkspace::gs_pop_get_children(gs1, y = "live")

# Subset count matrix to only populations of interest
keep_idx <- rownames(gs_counts) %in% subpops
clust_counts <- gs_counts[rownames(gs_counts) %in% subpops, ]

# Ensure columns are in order that matches design matrix
samples <- as.character(sample_info$filename)
clust_counts <- clust_counts[, samples]
rownames(clust_counts) <- sub("^.*/", "", rownames(clust_counts))
clust_counts
#>             sample_1.fcs sample_2.fcs sample_3.fcs sample_4.fcs sample_5.fcs
#> 26                     9            4           13            1            7
#> 27                   166          234          124           41           74
#> 28                    21           50           35           22           45
#> 29                   105          364          109          254           76
#> 30                    10           19            8           10           25
#> 31                    51          221           53          124           47
#> 32                    10           13           10           14           10
#> 33                    24           19           24           29           20
#> 34                     8           31           10           16            8
#> 35                    10           42           12           36           15
#> B cells               36           26          430           49           73
#> CD4 T cells          193          211          365          335          577
#> CD8 T cells           69           53          537         1122          805
#> Monocytes            368          812          912          291          786
#> NK T cells           127           35           74           40           23
#> NK cells              57           92          197           60           59
#> Undefined            115          106          372           58          100
#> cDC                   63            8            8            1            3
#> gdT cells              4            1           92            7            4
#> pDC                  123           43          304            6           18
#>             sample_6.fcs
#> 26                     8
#> 27                   117
#> 28                    57
#> 29                   129
#> 30                    33
#> 31                    71
#> 32                    10
#> 33                    27
#> 34                    13
#> 35                    28
#> B cells               61
#> CD4 T cells          336
#> CD8 T cells         1052
#> Monocytes            377
#> NK T cells            20
#> NK cells              25
#> Undefined            222
#> cDC                    1
#> gdT cells              7
#> pDC                   14
```

Given any matrix of counts,
[`plotSampleProportions()`](https://00berst33.github.io/flowFun/reference/plotSampleProportions.md)
visualizes cell type proportions across samples.

``` r

plotSampleProportions(clust_counts)
```

![](workflowVignetteUpdated_files/figure-html/visualize-counts-1.png)

#### Differential Abundance

Once these matrices have been acquired, to perform differential
abundance analysis, the user only needs to pass them to
[`doDAAnalysis()`](https://00berst33.github.io/flowFun/reference/doDAAnalysis.md),
and specify a normalization method, if any. Results are returned as a
list, where each entry is a table corresponding to a comparison. The
tables are named after the comparison they correspond to, and contain
results of likelihood ratio tests performed by edgeR, with each test
being ranked by its adjusted p-value.

``` r

# Perform differential abundance analysis
da_results <- doDAAnalysis(design = design, 
                           counts = counts, 
                           contrasts = contrasts,
                           sample_df = sample_info, 
                           norm_method = "TMM")
```

The resulting individual tables may be accessed with either double
brackets `[[]]` or `$`, e.g. `da_results[[1]]` or
`da_results$MIBC_NAC__vs__MIBC_No.NAC`. In this example, two comparisons
were tested, so the resulting list contains two tables:

|             |   logFC    |  logCPM  |     LR     |  PValue   |    FDR    |
|:------------|:----------:|:--------:|:----------:|:---------:|:---------:|
| cDC         | -3.7651269 | 12.67172 | 19.6367982 | 0.0000094 | 0.0001873 |
| CD8 T cells | 3.9582888  | 18.03725 | 17.0039740 | 0.0000373 | 0.0003730 |
| gdT cells   | 2.3350335  | 12.60229 | 5.2020767  | 0.0225599 | 0.1503995 |
| B cells     | 1.7467617  | 15.12860 | 4.1817193  | 0.0408622 | 0.2043111 |
| 27          | -1.3174968 | 15.59524 | 3.0938467  | 0.0785888 | 0.2691608 |
| NK T cells  | -1.3282430 | 14.42161 | 3.0497770  | 0.0807482 | 0.2691608 |
| pDC         | -1.1274898 | 14.71813 | 2.1883188  | 0.1390610 | 0.3973172 |
| CD4 T cells | 1.0344152  | 17.11514 | 1.6368888  | 0.2007529 | 0.5018822 |
| 31          | -0.7152660 | 15.30686 | 0.8858492  | 0.3466046 | 0.7702325 |
| 29          | -0.5955501 | 16.17769 | 0.6171876  | 0.4320943 | 0.8054287 |
| 34          | -0.6230041 | 12.69379 | 0.5885361  | 0.4429858 | 0.8054287 |
| Undefined   | 0.4693132  | 15.81039 | 0.3528109  | 0.5525262 | 0.9208770 |
| 33          | 0.2405540  | 13.37792 | 0.0870335  | 0.7679829 | 0.9984792 |
| 30          | 0.2046768  | 12.96451 | 0.0583372  | 0.8091436 | 0.9984792 |
| NK cells    | -0.1171831 | 14.88648 | 0.0227515  | 0.8801051 | 0.9984792 |
| 26          | -0.1109586 | 11.70795 | 0.0154018  | 0.9012329 | 0.9984792 |
| 28          | 0.0860343  | 13.97126 | 0.0115634  | 0.9143658 | 0.9984792 |
| Monocytes   | -0.0740140 | 17.80082 | 0.0093045  | 0.9231553 | 0.9984792 |
| 32          | -0.0047519 | 12.39884 | 0.0000361  | 0.9952036 | 0.9984792 |
| 35          | -0.0013819 | 13.41683 | 0.0000036  | 0.9984792 | 0.9984792 |

Ctrl vs. MIBC {.table}

|             |   logFC    |  logCPM  |    LR     |  PValue   |    FDR    |
|:------------|:----------:|:--------:|:---------:|:---------:|:---------:|
| pDC         | -2.3807569 | 14.71813 | 6.1353503 | 0.0132505 | 0.1584515 |
| gdT cells   | -2.4435551 | 12.60229 | 5.8199770 | 0.0158452 | 0.1584515 |
| 30          | 1.7488640  | 12.96451 | 3.1553614 | 0.0756779 | 0.5045195 |
| B cells     | -1.1126297 | 15.12860 | 1.5035665 | 0.2201234 | 0.6600647 |
| NK T cells  | -1.0937646 | 14.42161 | 1.3803578 | 0.2400402 | 0.6600647 |
| NK cells    | -1.0671431 | 14.88648 | 1.3648476 | 0.2426991 | 0.6600647 |
| 28          | 1.0721412  | 13.97126 | 1.3358235 | 0.2477719 | 0.6600647 |
| 29          | -1.0094303 | 16.17769 | 1.2475242 | 0.2640259 | 0.6600647 |
| 31          | -0.7542249 | 15.30686 | 0.6911796 | 0.4057634 | 0.8937357 |
| 27          | 0.6190904  | 15.59524 | 0.4697341 | 0.4931095 | 0.8937357 |
| 26          | 0.7049993  | 11.70795 | 0.4675487 | 0.4941170 | 0.8937357 |
| CD4 T cells | 0.5519135  | 17.11514 | 0.3825504 | 0.5362414 | 0.8937357 |
| Monocytes   | 0.4786227  | 17.80082 | 0.2887697 | 0.5910100 | 0.9092462 |
| cDC         | -0.5564371 | 12.67172 | 0.2205538 | 0.6386183 | 0.9123119 |
| 35          | -0.3650014 | 13.41683 | 0.1518197 | 0.6968024 | 0.9202165 |
| 34          | -0.3308023 | 12.69379 | 0.1135186 | 0.7361732 | 0.9202165 |
| 32          | -0.2124137 | 12.39884 | 0.0461951 | 0.8298216 | 0.9282212 |
| Undefined   | 0.1580130  | 15.81039 | 0.0313012 | 0.8595700 | 0.9282212 |
| 33          | -0.1387510 | 13.37792 | 0.0221042 | 0.8818101 | 0.9282212 |
| CD8 T cells | 0.0213681  | 18.03725 | 0.0005791 | 0.9808014 | 0.9808014 |

NAC vs. No NAC {.table}

#### Differential Expression

To test for differential expression, the design and contrast matrices
must again be specified, in addition to the channels/markers to be
tested. Recall that in most cases, these markers should NOT be the same
as those used for clustering.

Once the required information is prepared, performing a basic
differential expression analysis requires little input from the user, as
all necessary operations are wrapped up within the function
[`doDEAnalysis()`](https://00berst33.github.io/flowFun/reference/doDEAnalysis.md).
This function accepts either a `data.table` or `GatingSet`.

When calling the function, users have the option to give a count matrix,
which will be used to assign weights to each cluster-marker pair before
testing. We show two examples of this function below, called on our
`GatingSet` complete with all the nodes we defined during clustering.

``` r

# Define channels of interest
cols_to_test <- c("FITC-A", "BV711-A")

# Set subpopulations of interest to variable 
subpops <- gs_pop_get_children(gs1, 
                               "live", # this retrieves all clusters that are  
                                       # children of the "live" node
                               path = "auto") # this argument concatenates population names
# Remove `Undefined` from populations of interest
subpops <- subpops[!subpops %in% c("Undefined")]

# Perform differential expression analysis
de_res <- doDEAnalysis(gs1, 
                       cols_to_test = cols_to_test, 
                       design = design, 
                       contrasts = contrasts, # subset contrasts matrix to column of interest
                       counts = counts, # optional parameter, used to calculate weights when given
                       subpopulations = subpops,
                       print_table = TRUE) # TRUE by default
```

The output is a list, whose length corresponds to the number of
comparisons defined by `contrasts`. Each element is a formatted table
with features, average expression, and p-values (if desired, the raw
output from
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html) may be
returned instead of tables by setting `print_table=FALSE`). These tables
are printed neatly below.

|  Ctrl_X   | MIBC_NAC  | MIBC_No.NAC |
|:---------:|:---------:|:-----------:|
| 2.7752128 | 2.7036740 |  2.7764469  |
| 1.0948151 | 1.0876747 |  1.0940220  |
| 2.6113283 | 2.4868329 |  2.5938907  |
| 1.0536996 | 1.1585773 |  1.1509620  |
| 2.7507389 | 2.7325526 |  2.8957611  |
| 0.9946193 | 0.9294050 |  1.0108297  |
| 2.8915848 | 2.7020617 |  2.7733731  |
| 1.4614352 | 1.5826974 |  1.5323702  |
| 2.6627426 | 2.5957556 |  2.5334974  |
| 1.4809217 | 1.3882520 |  1.1975221  |
| 2.7326816 | 2.5473377 |  2.3569395  |
| 1.3079426 | 1.3258470 |  1.0696589  |
| 2.9515291 | 2.8736278 |  2.4367049  |
| 1.7174328 | 1.3779621 |  1.2203727  |
| 2.7698251 | 2.7476427 |  2.8978313  |
| 0.7539766 | 0.7863174 |  0.7755172  |
| 3.2217065 | 3.2008727 |  3.1705052  |
| 2.0523656 | 1.6696045 |  1.2155066  |
| 3.0548319 | 3.0679278 |  3.0620334  |
| 1.8890108 | 1.8791540 |  1.8417370  |
| 3.1867400 | 3.1468184 |  3.1774635  |
| 2.2133555 | 1.9726037 |  2.1364879  |
| 2.7537351 | 2.7578061 |  2.8794287  |
| 0.9923584 | 1.0605324 |  1.0993760  |
| 2.7567234 | 2.4945123 |  2.7622536  |
| 1.7826523 | 1.4854272 |  1.0730557  |
| 3.1837451 | 3.1583203 |  3.1699445  |
| 2.3754419 | 2.2301533 |  2.3121723  |
| 2.5886083 | 2.4623907 |  2.5888142  |
| 2.7939553 | 2.3351218 |  2.3359634  |
| 2.8603506 | 2.6947581 |  2.7318391  |
| 2.0184685 | 1.7971035 |  1.7432754  |
| 2.9831529 | 2.7625465 |  2.9035543  |
| 2.6804596 | 2.5042479 |  2.5537776  |
| 2.8671712 | 2.7111781 |  2.8241923  |
| 2.4105255 | 2.1971313 |  2.3378892  |
| 2.7847634 | 2.7573784 |  2.9658266  |
| 1.1125609 | 1.4058297 |  1.2622483  |

Ctrl vs. MIBC {.table}

|  Ctrl_X   | MIBC_NAC  | MIBC_No.NAC |
|:---------:|:---------:|:-----------:|
| 0.0117234 | 0.0117234 |  0.0117234  |
| 0.0117234 | 0.0117234 |  0.0117234  |
| 0.0272166 | 0.0272166 |  0.0272166  |
| 0.0272166 | 0.0272166 |  0.0272166  |
| 0.0395904 | 0.0395904 |  0.0395904  |
| 0.0395904 | 0.0395904 |  0.0395904  |
| 0.0118745 | 0.0118745 |  0.0118745  |
| 0.0118745 | 0.0118745 |  0.0118745  |
| 0.0313728 | 0.0313728 |  0.0313728  |
| 0.0313728 | 0.0313728 |  0.0313728  |
| 0.0319438 | 0.0319438 |  0.0319438  |
| 0.0319438 | 0.0319438 |  0.0319438  |
| 0.0771517 | 0.0771517 |  0.0771517  |
| 0.0771517 | 0.0771517 |  0.0771517  |
| 0.1091089 | 0.1091089 |  0.1091089  |
| 0.1091089 | 0.1091089 |  0.1091089  |
| 0.0257172 | 0.0257172 |  0.0257172  |
| 0.0257172 | 0.0257172 |  0.0257172  |
| 0.0466252 | 0.0466252 |  0.0466252  |
| 0.0466252 | 0.0466252 |  0.0466252  |
| 0.0219581 | 0.0219581 |  0.0219581  |
| 0.0219581 | 0.0219581 |  0.0219581  |
| 0.0157446 | 0.0157446 |  0.0157446  |
| 0.0157446 | 0.0157446 |  0.0157446  |
| 0.0690066 | 0.0690066 |  0.0690066  |
| 0.0690066 | 0.0690066 |  0.0690066  |
| 0.0296957 | 0.0296957 |  0.0296957  |
| 0.0296957 | 0.0296957 |  0.0296957  |
| 0.0863868 | 0.0863868 |  0.0863868  |
| 0.0863868 | 0.0863868 |  0.0863868  |
| 0.0591312 | 0.0591312 |  0.0591312  |
| 0.0591312 | 0.0591312 |  0.0591312  |
| 0.0762493 | 0.0762493 |  0.0762493  |
| 0.0762493 | 0.0762493 |  0.0762493  |
| 0.0591312 | 0.0591312 |  0.0591312  |
| 0.0591312 | 0.0591312 |  0.0591312  |
| 0.0659380 | 0.0659380 |  0.0659380  |
| 0.0659380 | 0.0659380 |  0.0659380  |

NAC vs. No NAC {.table}

[`doDEAnalysis()`](https://00berst33.github.io/flowFun/reference/doDEAnalysis.md)
creates a table of cluster-sample MFIs internally, but some users may
wish to have this table for use outside of differential analysis. Two
functions are available for this purpose. If expression data is in a
`data.table` format, the function
[`getExprMatDE()`](https://00berst33.github.io/flowFun/reference/getExprMatDE.md)
may be used. If using a `GatingSet` object, the function
[`gs_makeMFIMatrix()`](https://00berst33.github.io/flowFun/reference/gs_makeMFIMatrix.md)
should be used instead. The result is a table with rows corresponding to
each cluster-marker pair, and columns for each sample.

Below we show an example of the output from
[`gs_makeMFIMatrix()`](https://00berst33.github.io/flowFun/reference/gs_makeMFIMatrix.md).
Note that `inverse = TRUE` is only applicable for `GatingSets` with
associated transformations;
[`getExprMatDE()`](https://00berst33.github.io/flowFun/reference/getExprMatDE.md)
is faster but unable to back-transform data.

``` r

# Print table
gs_makeMFIMatrix(gs1,
                 cols = cols_to_test,
                 subpopulations = subpops,
                 inverse = TRUE) %>%
  head()
#> # A tibble: 6 × 7
#>   feature       sample_1.fcs sample_2.fcs sample_3.fcs sample_4.fcs sample_5.fcs
#>   <chr>                <dbl>        <dbl>        <dbl>        <dbl>        <dbl>
#> 1 CD8 T cells.…       4922.        5280.         3817.       6861.        3146. 
#> 2 CD8 T cells.…        161.         150.          155.        155.         178. 
#> 3 B cells.FITC…       3776.        3320.         3660.       3171.        1858. 
#> 4 B cells.BV71…        165.          98.6         127.        256.         297. 
#> 5 NK T cells.F…       4238.        5499.         3609.      12470.        3514. 
#> 6 NK T cells.B…         94.3        101.          123.         90.3         62.6
#> # ℹ 1 more variable: sample_6.fcs <dbl>
```

Difference in marker expression between groups may be visualized with
the function
[`plotGroupMFIBars()`](https://00berst33.github.io/flowFun/reference/plotGroupMFIBars.md).
As above, if input is a `GatingSet` with transformations attached, the
data may be back-transformed to the linear scale before the calculation
of MFIs or statistical tests.

``` r

# Generate bar plot
plotGroupMFIBars(gs1,
                 col = "FITC-A",
                 sample_df = sample_info, 
                 comparison = comparisons[[1]],
                 populations = subpops,
                 inverse = TRUE)
```

![](workflowVignetteUpdated_files/figure-html/plot-mfis-1.png)

Finally, grouped dimension reduction plots may be generated, colored
either by density or marker expression, then compared against their
parent plot if one exists.

``` r

# Make new parent umap
umap <- plotUMAP(fsom_dt, num_cells = 7500, seed = 42)

# Color by density
plotGroupUMAPs(fsom_dt, sample_info, comparisons[[1]], umap = umap)

# May also color by marker expression
plotGroupUMAPs(fsom_dt, sample_info, comparisons[[1]], color_by = "BV711-A", umap = umap)
```

![](images/workflow_final/group_umap1.png)![](images/workflow_final/group_umap2.png)![](images/workflow_final/parent_umap.png)
