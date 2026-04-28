# Perform clustering on a group of preprocessed flow cytometry data files. The
# clustering is done with the package FlowSOM, which works by constructing
# a self-organizing map (SOM), a type of artificial neural network. This SOM is
# created in a way such that the closer two nodes are to each other, the more
# similar they are to each other. Through the process of constructing this SOM,
# each cell in the data is mapped to whichever node (or cluster) in the SOM that
# represents it best, resulting in the final clustering. Typically, depending on
# the desired resolution, a very high number of clusters is used - anywhere from
# 50 - 200 may be reasonable. The clustering is then visualized as a minimum
# spanning tree (MST).
#
# After clustering is complete, the next step is to identify cell type
# populations. Because we create such a large number of nodes in the first step,
# it is necessary to cluster them further so that we may interpret our
# results in a biologically meaningful way. It is possible to use FlowSOM to
# cluster the nodes, essentially producing clusters of clusters, which we
# call "metaclusters". For example, if we expected to see CD4 T cells,
# CD8 T cells, NK T cells, and gdT cells in our data, we might tell FlowSOM that
# we want to first cluster our cells into 100 nodes, then cluster those 100 nodes
# into 4 metaclusters. However, FlowSOM lacks any biological knowledge about
# our data, and does not always reliably distinguish cell types. Because of this,
# this workflow uses the strategy of overclustering, and requires the user to
# manually merge metaclusters by examining dimension reduction plots and heatmaps.
#
# Finally, this script allows the user to perform backgating and reclustering on
# any identified cell types of particular interest. This clustering may be done
# using either a selected number of principal components obtained from PCA,
# or a new subset of markers defined by the user.

# Load packages.
library(data.table)
library(flowFun)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(ggplot2)

# Set a seed for reproducibility
set.seed(84)

# Set the directory you would like to store analysis results in.
work_dir <- file.path("C:/Users/00ber/OneDrive/Desktop/VPC")

#####
## Load in raw data for first time
# The name of the directory containing the .fcs files to analyze
data_dir <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw" ###

# Get all filenames and create GatingSet
files <- list.files(data_dir, full.names = TRUE)
cs <- flowWorkspace::load_cytoset_from_fcs(files,
                                           which.lines = 5000) # if NULL, all cells are read in
                                                                # if an integer, a random sample of given size is read in
gs <- flowWorkspace::GatingSet(cs)

###
## OR, Load in existing data
# note: changes made to object created from load_cytoset() do not affect
#   original data unless `backend_readonly=FALSE`
# Load previously created GatingSet or cytoset, if needed
gs <- flowWorkspace::load_gs(file.path(work_dir, "template_gs"))
# cs <- flowWorkspace::load_cytoset(file.path(work_dir, "template_cs"))
###

# Make deep clone of GatingSet, to avoid making accidental changes to underlying data
gs1 <- flowWorkspace::gs_clone(gs)

# Compensate
comp_mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/morgans_comp_matrix.csv",
                     check.names = FALSE) # note check.names

# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.
# Check
colnames(comp_mat)
colnames(gs1)

# Print column names in GatingSet not found in column names of compensation matrix
colnames(gs1)[!(colnames(gs1) %in% colnames(comp_mat))]
# Prepare compensation matrix to be passed to `compensate`
#   Columns are matched to those in the GatingSet based on the regexpr given to pattern.
#   Default is `" >.*"`. Setting `pattern = ""` tells the function that column names
#   do not need to be matched.
comp_mat <- prepareCompensationMatrix(comp_mat, gs1, pattern = "")

# Compensate data
compensate(gs1, comp_mat)

# Define custom transformation, if desired
# logicle transformation, standard for flow cytometry:
log_trans <- flowCore::estimateLogicle(gs1[[1]],
                                   channels = colnames(comp_mat))
# Or inverse hyperbolic sin transform, standard for mass cytometry:
asinh_trans <- flowWorkspace::asinhtGml2_trans(equal.space = TRUE)
asinh_trans <- flowWorkspace::transformerList(colnames(comp_mat), asinh_trans)

# getting inverse
# asinh_inv <- flowCore::transformList(names(asinh_trans), lapply(asinh_trans, `[[`, "inverse"))


# Apply transformation
flowWorkspace::transform(gs1, log_trans)

# Specify stain for dead cells, if you have one
ld_stain <- "BUV496-A"

# Prepare gatingTemplate for preprocessing
gt_table <- generateGatingTable(gs1,
                                collapse_data = TRUE, # if TRUE, gates are drawn on collapsed data and replicated across all samples
                                ld_stain = "BUV496-A") # remove this argument if your data has no L/D stain

# Check resulting table
gt_table

# Make changes if desired
# For example, to change the alias of the population "live", specify the corresponding row
#   and column name in brackets:
gt_table[4, "alias"] <- "live_cells"

# Check results
gt_table

# Turn table into gatingTemplate once satisfied
gt <- openCyto::gatingTemplate(gt_table)

# Visualize gating template
openCyto::plot(gt)

# Apply gating scheme
gt_gating(gt, gs1)

## Make graphs to check results of preprocessing:

## Check one gate for each sample
# nonDebris gate
plotAllSamples(gs1, "FSC-A", "SSC-A", "nonMargins", "nonDebris")

# singlets gate
plotAllSamples(gs1, "FSC-A", "FSC-H", "nonDebris", "singlets")

# live cell gate
plotAllSamples(gs1, !!enquo(ld_stain), "FSC-A", "singlets", "live_cells")


# Adjust gates if necessary, either by editing and reapplying the gating template,
#   or manually redrawing them.

### Edit arguments of gating template
gt_table[2, "gating_args"] <- "gate_range=c(0,80000)"

# Remove gate(s) that will be redrawn
# Here we delete "nonDebris", which also deletes its children gates (singlets, live_cells)
flowWorkspace::gs_pop_remove(gs1, "nonDebris")

# Generate new gating template
gt <- openCyto::gatingTemplate(gt_table)

# Apply to data
openCyto::gt_gating(gt, gs1)

### Redraw gates
library(CytoExploreR)

# Set desired name to save gating template under
gt_name <- "gatingTemplate.csv"

# Write gatingTemplate to CSV file
write.csv(gt_table, file.path(work_dir, gt_name))

# Manually redraw desired gate
CytoExploreR::cyto_gate_edit(gs1,
                             parent = "singlets", # parent population
                             alias = "live_cells", # name of the gate to edit
                             channels = c("BUV496-A", "FSC-A"), # channels to gate on
                             type = "polygon", # type of gate to draw
                             gatingTemplate = file.path(work_dir, gt_name))


# Once you have results, save GatingSet
flowWorkspace::save_gs(gs1, path = file.path(work_dir, "template_gs"))


##########
# CytoML #
##########

# Open FlowJo workspace in R from .xml file
flowjo_file <- "path/to/flowjo.xml"
ws <- CytoML::open_flowjo_xml(file)

# Make GatingSet from FlowJo workspace
gs <- CytoML::flowjo_to_gatingset(ws,
                                  path = "path/to/fcs_files") # or cytoset = ...

# Get first sample
gh <- gs[[1]]

# Plot gates of first sample
ggcyto::autoplot(gh)


# Must download docker image
# To save current GatingSet as a FlowJo workspace;
CytoML::gatingset_to_flowjo(gs, "path/to/saved_ws.wsp")


##############
# Clustering #
##############

# Get flowSet from GatingSet
ex_fs <- flowWorkspace::gs_pop_get_data(gs1, "live_cells")
ex_fs <- flowWorkspace::cytoset_to_flowSet(ex_fs)

# Perform clustering
# If you plan to apply controls to your data and calculate delta MFIs, it is
#   highly recommended that you save the FlowSOM object at this stage by setting
#   `fsom_file` equal to a filepath.
fsom_dt <- flowSOMWrapper(ex_fs,
                          xdim = 10,
                          ydim = 10,
                          cols_to_cluster = c(10, 12:14, 16, 18:23, 25:32, 34), # Define markers/columns to use for clustering
                          num_clus = 23,
                          seed = 42,
                          fsom_file = "fsom.rds") # if you intend to use controls, it is helpful to save the clustering

### Iteratively merge clusters until all cell types identified
# Generate heatmap
plotMetaclusterMFIs(fsom_dt)
# Generate UMAP
plotUMAP(fsom_dt, num_cells = 2500, seed = 42)
# Generate 2D scatterplot
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV563-A"),
                     clusters = c(4, 42, 76),
                     metaclusters = NULL)

# Merge metaclusters
fsom_dt <- editTableMetaclusters(fsom_dt,
                                 new_labels = c("9" = "1",
                                                "2" = "1",
                                                "8" = "4",
                                                "23" = "15",
                                                "21" = "15",
                                                "20" = "15",
                                                "22" = "15",
                                                "7" = "6"))

# Generate annotated heatmap when final clustering reached
annotateMFIHeatmap(fsom_dt, cols_to_cluster)


### Final clusters may be added to GatingSet as a boolean gate
# See all populations in GatingSet, choose gate whose cells were used for clustering
flowWorkspace::gs_get_pop_paths(gs1)

# Add gates for each cluster to GatingSet
addClustersToGatingSet(fsom_dt,
                       gs1,
                       parent_gate = "live_cells",
                       fsom_file = NULL) # by default, this argument is NULL and the
                                         # function checks the Gating

# Visualize new gating template
openCyto::plot(gs1)

#############
# DA and DE #
#############

# Read in .csv file containing info about samples (must include a column for filename and one for sample name)
file <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv"

# Define list of comparisons we would like to make between groups
comparisons <- list(
  ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

# Prepare metadata for further analysis
sample_info <- prepareSampleInfo(file,
                                 name_col = "Sample.Name",
                                 filename_col = "File.Name",
                                 comparisons = comparisons)

# Generate design matrix
design <- makeDesignMatrix(sample_info)

# Generate contrasts matrix
contrasts <- makeContrastsMatrix(sample_info, comparisons)

# Generate matrix of sample/metacluster cell counts
counts <- makeCountMatrix(fsom_dt,
                          #meta_names = meta_of_interest,
                          min_cells = 3,
                          min_samples = 4)

#
### Perform differential abundance analysis
#
da_results <- doDAAnalysis(design = design,
                           counts = counts,
                           contrasts = contrasts,
                           sample_df = sample_info,
                           norm_method = "TMM")

#
### Perform differential expression analysis
#

# Set markers of interest
marker_cols <- c("BV711-A", "FITC-A")

# Get clustered populations
subpops <- gs_pop_get_children(gs1, "live_cells")

# Perform differential expression analysis for given markers
de_res <- doDEAnalysis(gs1,
                       cols_to_test = marker_cols,
                       design = design,
                       contrasts = contrasts,
                       subpopulations = subpops,
                       inverse = FALSE)

# View results
limma::topTable(de_res)

# Plots
# Get table with MFIs where rows are sample and columns are metaclusters
plot_mat <- getSampleMetaclusterMFIs(fsom_dt, "BV711-A", sample_info)

# Generate bar plot
plotGroupMFIBars(plot_mat,
                 sample_df = sample_info,
                 comparison = comparisons[[1]])

#### Optionally, apply 1D boundary gates after clustering
# Plot channel marker densities by sample/metacluster
plot1DMarkerDensities(gs1,
                      channel = "FITC-A",
                      population = "live",
                      facet_by = "subpopulations", # may also facet by "samples"
                      inverse = FALSE)

# Add gate to GatingSet with openCyto
# Add boundary gate on channel; note arguments here are the same as column names for
#   the earlier gatingTemplate
openCyto::gs_add_gating_method(gs1,
                               alias = "+FITC-A",
                               parent = "9",         # add +FITC-A gate to metacluster 9
                               dims = "FITC-A",
                               gating_method = "gate_mindensity",
                               collapseDataForGating = TRUE,
                               groupBy = length(gs1))

####
## Optionally, apply controls to find delta MFIs

# Load in control FCS files
ctrl_dir <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/FMO"

# Read in
ctrl <- flowWorkspace::load_cytoset_from_fcs(list.files(ctrl_dir, full.names = TRUE), which.lines = 20000)
# Apply transformations, compensations and gates to control gs
ctrl_gs2 <- flowWorkspace::gh_apply_to_cs(gs1[[1]], ctrl, compensation_source = "template") # make sure to exclude boolean

ctrl_fs <- flowWorkspace::cytoset_to_flowSet(ctrl)

# Apply clustering to controls
# This may also be done manually using `FlowSOM::NewData`
fsom_projected <- clusterControls(ctrl_gs, gs1, "live")

# Edit clusters if desired
# ...

# Get data.table for controls, and add clusters to corresponding GatingSet
ctrl_dt <- flowSOMToTable(fsom_projected)
addClustersToGatingSet(ctrl_dt, ctrl_gs, "live")

# Read in table of sample info, if you have one; see `?addMetadataToGatingSet` for
# details on how this table should be defined
# Should have column named 'filename'
sample_info <- read.csv("path/to/sample/info/csv")
# Add metadata to GatingSet
addMetadataToGatingSet(gs1, sample_info)
# Check results
pData(gs1)
# this also allows us to use functions in ggcyto package to facet by any group defined
# in this metadata

# Define markers to test
cols_to_test <- c("PHA-L", "IL10R")

# Get delta MFIs
delta_mfis <- gs_makeDeltaMFIs(gs1, ctrl_gs, subpopulations = subpops, metadata_col = "FMO")

# Do DE testing and make plots as you would for typical DE testing
