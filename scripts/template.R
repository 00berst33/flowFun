# This template includes all the steps to run the completed flowFun pipeline, 
# starting with FCS files through to statistical analysis and plotting of data. 

# Vignettes and other helpful documentation, including a "Getting Started with
# flowFun Guide", may be viewed on the package's main page: 
# https://00berst33.github.io/flowFun

# The vignette can also be opened in RStudio by running:
vignette(package = "flowFun", "workflowVignetteUpdated")

# Load packages.
library(data.table)
library(flowFun)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(ggplot2)
library(CytoExploreR) # I added this here but we could remove it and just keep it below if you think that's better


# Set a seed for reproducibility
set.seed(84)

# Set the directory you would like to store analysis results in.
work_dir <- file.path("/Users/morganroberts/Library/CloudStorage/OneDrive-VancouverProstateCentre/Data/Bioinformatics/High Parameter Flow Data/Testing FlowFun/April 2026/26-04-16_Analysis")

##################
# PRE-PROCESSING #
##################

### LOAD FCS FILES ###

# Note: Instead of FCS files, you can load a previously created GatingSet or 
# cytoset, if needed. See the "Getting Started with flowFun Guide" for instructions
# on how to do that.

# Indicate the name of the directory containing the .fcs files to analyze
data_dir <- "/Users/morganroberts/Library/CloudStorage/OneDrive-VancouverProstateCentre/Data/FACS Data/24-12-18 SL012 Aged Mice Blood Spl Bladder/SL012 Spleen Samples"

# Get all filenames and create a GatingSet
files <- list.files(data_dir,
                    pattern = "\\.fcs$",
                    full.names = TRUE,
                    ignore.case = TRUE)
cs <- flowWorkspace::load_cytoset_from_fcs(files,
                                           which.lines = 50000) # if NULL, all cells are read in
                                                                # if an integer, a random sample of given size is read in
gs <- flowWorkspace::GatingSet(cs)

# Make deep clone of GatingSet, to avoid making accidental changes to underlying data
gs1 <- flowWorkspace::gs_clone(gs)

### APPLY COMPENSTATION ####

# Load the compensation matrix
comp_mat <- read.csv("/Users/morganroberts/Library/CloudStorage/OneDrive-VancouverProstateCentre/Data/FACS Data/24-12-18 SL012 Aged Mice Blood Spl Bladder/SL012_Comp_Matrix_noQC_fixed_names.csv",
                     check.names = FALSE) # note check.names

# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.

# Check column names
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
# NOTE: If you plan to apply controls to your data and calculate delta MFIs, it is
#   highly recommended that you save the FlowSOM object at this stage by setting
#   `fsom_file` equal to a filepath. If you don't do this, at the very least
#   make sure to specify a seed and make note of the parameters you use in this
#   call, so that this particular clustering may be recreated later if needed.
fsom_dt <- flowSOMWrapper(ex_fs,
                          cols_to_cluster = c(10, 12:14, 16, 18:23, 25:32, 34), # Define markers/columns to use for clustering
                          num_clus = 23,
                          xdim = 10,
                          ydim = 10,
                          seed = 42,
                          fsom_file = "fsom.rds") # the resulting FlowSOM object will be saved to disk under this name

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
                                                "2" = "1", # MC (metacluster) 9 and 2 are merged into MC 1
                                                "8" = "4", # MC 8 is merged in MC 4
                                                "23" = "15",
                                                "21" = "15",
                                                "20" = "15",
                                                "22" = "15", # MC 23, 21, 20, and 22 are merged into MC 15
                                                "7" = "6")) # MC 7 is merged into MC 6

# This function can also be used to renamed metaclusters, once you are confident
# in a more biologically meaningful cell type
fsom_dt <- editTableMetaclusters(fsom_dt,
                                 new_labels = c("1" = "Monocytes",
                                                "4" = "NK cells",
                                                "15" = "CD8 T cells",
                                                "6" = "Undefined"))


# Generate annotated heatmap when final clustering reached
annotateMFIHeatmap(fsom_dt, cols_to_cluster)


### Final clusters may be added to GatingSet as a boolean gate
# See all populations in GatingSet, choose gate whose cells were used for clustering
flowWorkspace::gs_get_pop_paths(gs1)

# Add gates for each cluster to GatingSet
addClustersToGatingSet(fsom_dt,
                       gs1,
                       parent_gate = "live_cells",
                       fsom_file = NULL) # by default, this argument is NULL

###
## The following notes are relevant only if you intend to apply controls later:

# If you specified the parameter `fsom_file` when calling `flowSOMWrapper()`
# earlier, it is fine to leave `fsom_file` `NULL` here. `addClustersToGatingSet()`
# will find the file by checking the attributes of its input, in this case `fsom_dt`.
# The original FlowSOM object's metaclusters will be edited and renamed according to
# any edits you made, and its filename will be associated with your GatingSet,
# making later clustering of controls straightforward.

# If you did not specify `fsom_file` earlier, or wish to use a different file,
# the correct file should be specified in this call to `addClustersToGatingSet()`.
# The RDS file given must already exist.

# If the RDS file is successfully accessed, its name will be added to the metadata of the
# GatingSet. More specifically, it will be added to the data frame viewed with `pData(gs1)`,
# under a column name corresponding to the population/node the clustering was done on
# (in this case "live_cells").
###

# Visualize new gating template
openCyto::plot(gs1)

#############
# DA and DE #
#############

# Read in a .csv file containing info about samples
# IMPORTANT: Each row must correspond to a sample. The file must include a
# column for filename, and one for sample name. Any information about experimental
# group and/or corresponding control files that are relevant to differential
# analysis should also be included in the table; see the workflow vignette for
# more details.
file <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv"

# Define list of comparisons we would like to make between groups.
# Type `?prepareSampleInfo` in the console and scroll down to the examples to
# get more information on how `comparisons` should be defined.
comparisons <- list(
  ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

# Prepare metadata for further analysis. Returns a table read in from the given
# CSV file, with an added column called `group`, which is used to ensure that
# the correct samples are used in each test. Additionally, the column names given
# to `name_col` and `filename_col` are renamed to "sample.name" and "filename",
# respectively.
sample_info <- prepareSampleInfo(file,
                                 name_col = "Sample.Name", # name of column containing sample names in your file
                                 filename_col = "File.Name", # name of column containing file names in your file
                                 comparisons = comparisons)

# Generate design matrix
design <- makeDesignMatrix(sample_info)

# Generate contrasts matrix
contrasts <- makeContrastsMatrix(sample_info, comparisons)

# Select populations to find counts for, if you would like to exclude any cell types found earlier
meta_of_interest <- c("Monocytes", "NK cells", "CD8 T cells")

# Generate matrix of sample/metacluster cell counts
counts <- makeCountMatrix(fsom_dt,
                          populations = meta_of_interest, # if `fsom_dt` is a data frame, the default is all populations
                          min_cells = 3,
                          min_samples = 4)

# If desired, this matrix may be exported as a .csv file (the same goes for any other matrix or data.frame in R)
write.csv(counts, "table_filename.rds")

#
### Perform differential abundance analysis
# doDAAnalysis() returns a list of data.frames, where each element corresponds to
# a column in `contrasts` (i.e. a comparison)
da_results <- doDAAnalysis(design = design,
                           counts = counts,
                           contrasts = contrasts,
                           sample_df = sample_info,
                           norm_method = "TMM")

# Results may be viewed by either simply typing `da_results` in the console,
# or selecting a particular comparison in this list to view. The example `comparisons`
# object defined above had two elements, so we could view results as follows:
da_results
# or
da_results[[1]] # display only first comparison, MIBC vs. Ctrl
da_results[[2]] # display only second comparison, NAC vs. No NAC

# You will see five columns in the resulting tables: logFC, logCPM, LR, PValue, FDR
# i.e., log-fold change, log counts per million, likelihood ratio, p-value, false discovery rate (essentially adjusted p-value)
# The columns PValue and FDR will likely be most notable, as the first three columns are more relevant to RNAseq data.

#
### Perform differential expression analysis
#

# Set markers of interest
# IMPORTANT: These should NOT be markers used for clustering
marker_cols <- c("BV711-A", "FITC-A")

# Get clustered populations
subpops <- gs_pop_get_children(gs1, "live_cells", path = "auto")

# Perform differential expression analysis for given markers. This function accepts
# either a GatingSet or data.table.
de_res <- doDEAnalysis(gs1,
                       cols_to_test = marker_cols,
                       design = design,
                       contrasts = contrasts,
                       subpopulations = subpops, # only necessary to specify when input is a GatingSet
                       inverse = TRUE) # when input is a GatingSet and inverse is TRUE,
                                       # data is back-transformed before testing

# View results
limma::topTable(de_res)

# Plots
# Generate MFI bar plot
plotGroupMFIBars(gs1,
                 col = "BV711-A", # name of channel to plot MFIs for
                 sample_df = sample_info,
                 comparison = comparisons[[1]],
                 populations = subpops, # name of populations to plot
                 inverse = TRUE, # only applicable if input is GatingSet
                 upper_lim = NULL)

# If you would like to get a table of sample/metacluster MFIs to use for your own plotting,
# use `getSampleMetaclusterMFIs()`.
plot_mat <- getSampleMetaclusterMFIs(gs1, "BV711-A", sample_df = sample_info,
                                     populations = subpops, inverse = TRUE)
write.csv(plot_mat, "tim3_mfis.csv")

#### Optionally, apply 1D boundary gates after clustering
# Plot channel marker densities by sample/metacluster
plot1DMarkerDensities(gs1,
                      channel = "FITC-A",
                      population = "live_cells",
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

# Specify directory containing the control FCS files you'd like to apply
ctrl_dir <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/FMO"

# Read FCS files into a cytoset
ctrl <- flowWorkspace::load_cytoset_from_fcs(list.files(ctrl_dir, full.names = TRUE), which.lines = 20000)
# Apply transformations, compensations and gates created above for main data (the `gs1` object)
# to this cytoset. Essentially, we are pre-processing the control data exactly
# as we did for our data earlier
ctrl_gs <- flowWorkspace::gh_apply_to_cs(gs1[[1]], ctrl, compensation_source = "template")

# Check that the gates are appropriate and don't need adjustment
# non-debris gate
plotAllSamples(ctrl_gs, "FSC-A", "SSC-A", "nonMargins", "nonDebris")
# singlets gate
plotAllSamples(ctrl_gs, "FSC-A", "FSC-H", "nonDebris", "singlets")


# Cluster control samples using same mapping as main data. `clusterControls()` returns a FlowSOM object.
# NOTE: This function assumes that the object passed to `primary_gs`, in this
# case `gs1`, has an associated FlowSOM object from calling `addClustersToGatingSet()` earlier.
# You can check for any associated FlowSOM RDS files with flowWorkspace::pData(gs1).
fsom_projected <- clusterControls(ctrl_gs, gs1, "live_cells")

# Get data.table for controls, and add clusters to corresponding GatingSet
ctrl_dt <- flowSOMToTable(fsom_projected)
addClustersToGatingSet(ctrl_dt, ctrl_gs, "live")

# Read in table of sample info, if you have one. This may be the result of
# `prepareSampleInfo()`, like the data frame created above for differential analysis.
# MUST have a column called "filename".
sample_info <- read.csv("path/to/sample/info/csv")

# Add metadata to GatingSet
addMetadataToGatingSet(gs1, sample_info)
# Check results, the columns in `sample_info` should be added to `pData(gs1)`
pData(gs1)
# (This also allows us to use functions in ggcyto package to facet by any group defined
# in this metadata.)

# Define markers to apply controls to
cols_to_test <- c("PHA-L", "IL10R")

# Get delta MFIs
delta_mfis <- gs_makeDeltaMFIs(gs1,
                               ctrl_gs,
                               subpopulations = subpops,
                               cols = cols_to_test,
                               metadata_col = "FMO") # the name of the column containing control filenames
                                                     # (must be present in pData(gs1))

# Do DE testing and make plots as you would for typical DE testing
