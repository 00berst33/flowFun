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
cs <- flowWorkspace::load_cytoset_from_fcs(files)
gs <- flowWorkspace::GatingSet(cs)

# Save
flowWorkspace::save_gs(gs, path = file.path(work_dir, "example_gs"))

## OR, Load in existing data
# Load previously created GatingSet, if needed
# gs <- flowWorkspace::load_gs(file.path(tmp, "example_gs"))

# Make deep clone of GatingSet
gs1 <- flowWorkspace::gs_clone(gs)

# Compensate
comp_mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/morgans_comp_matrix.csv",
                     check.names = FALSE) # note check.names

# Define custom transformation, if desired
# Example:
# trans <- flowCore::estimateLogicle(gs[[1]],
#                                    channels = colnames(comp_mat))

# Specify stain for dead cells, if you have one
ld_stain <- "BUV496-A"

# Preview preprocessing steps
previewPreprocessing.GatingSet(gs1,
                               sample_ind = 1,
                               compensation = comp_mat,
                               transformation = "logicle", # may be a custom `transformList`
                               ld_channel = ld_stain,
                               debris_args = list(),
                               singlet_args = list(),
                               live_args = list())

# Apply preprocessing steps
doPreprocessing.GatingSet(gs1,
                          compensation = comp_mat,
                          transformation = "logicle",
                          ld_channel = ld_stain,
                          debris_args = "gate_range = c(0, 75000)",             # see also `openCyto::gate_mindensity()`
                          singlet_args = list(),
                          live_args = "gate_range = c(0, 2)")                   # see also `flowStats::gate_singlet()`

## Make graphs to check results of preprocessing:
# Non-debris gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `SSC-A`), subset = "nonDebris") + #add subset attribute before here?
  geom_hex(bins = 200) +
  geom_gate("nonDebris") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Singlet gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `FSC-H`), subset = "nonDebris") +
  geom_hex(bins = 200) +
  geom_gate("singlets") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Live cell gate
ggcyto::ggcyto(gs1, mapping = aes(x = !!enquo(ld_stain), y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 200) +
  geom_gate("live") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Once you have results, save and convert data to flowSet to perform clustering
flowWorkspace::save_gs(gs1, path = file.path(work_dir, "preprocessed_gs")) ### should cytoset be saved instead?
fs <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs1))

# After clustering; how to save results? write to cytoset? how to gate on subsets?

##############
# CLUSTERING #
##############

# Run FlowSOM algorithm. It typically only takes a few minutes to run, but note
# that runtime and memory usage are linearly related to the number of cells
# used in the clustering. So, if you decide to double the amount of cells used,
# the algorithm will take twice as long to run as it did before, and use twice
# as much memory.
#
# By default, 100 clusters are created, and 30 metaclusters
# are created. If you would like a higher resolution, you may use more clusters
# by editing the "xdim" and "ydim" parameters, which specify the dimensions of
# the SOM. You may also edit the number of metaclusters used by editing the
# "num_clus" variable below.

# Define the markers that you would like to use for clustering. Please make sure
# that you type the marker names the same as they appear in the .fcs files.
# IMPORTANT: These markers should not overlap with those that you would like to
# test for differential expression.

# May also be defined as a numeric vector of indices
cols_to_cluster <- c(12, 14:16, 18, 20:25, 27:34, 36)

## Perform clustering
# Returns FlowSOM object
fsom <- FlowSOM::FlowSOM(fs,
                         colsToUse = cols_to_cluster,
                         seed = 42)

# Returns tidytable
fsom_dt <- flowSOMWrapper(fs,
                          cols_to_cluster = cols_to_cluster,
                          num_clus = 30,
                          xdim = 10,
                          ydim = 10,
                          seed = 42) # ??? vs. seed above

# The created FlowSOM object is automatically saved to the current working
# directory. You may edit the file name with the parameter `fsom_dt_file`.

# Alternatively, if you have a previously created FlowSOM object that you would
# like to continue analysis of, you may read it in by uncommenting and running
# the line below with the filename of your FlowSOM object.
fsom_dt <- readRDS("fsom_A.rds")

# In addition to the star plots, dimension reduction plots and heatmaps are
# critical in determining the cell type of each metacluster, and which of said
# metaclusters should be merged. Functions to create these plots are defined below.
#
# Note - when deciding whether or not to merge metaclusters, your main reference
# is the heatmap of marker expression across metaclusters. The dendrogram clustering
# the metaclusters is of particular interest, and this hierarchy may be followed to
# perform the clustering. However, this dendrogram should not be followed blindly.
# First, the metaclusters that are determined to be similar by the clustering may
# not be of biological interest; some cell type markers may be of greater interest
# than others, and this is something that the clustering does not take into
# account. Second, dendrograms may be created with different linkage methods
# (single-linkage, average-linkage, centroid-linkage, etc.), and these different
# methods may result in slightly different hierarchies. This is not to say that
# the dendrogram is useless - only to inform the user of its shortcomings.
#
# Furthermore, the user should also consider cluster size when merging. If a
# metacluster is rather small, and differs from its neighboring clusters in
# irrelevant markers, then it is reasonable to merge it. However, if the
# metacluster is small but distinct in markers of interest, it may be left
# unmerged. Finally, if a metacluster does not represent any particular cell
# type of interest, it may be dropped from the analysis entirely. For example,
# depending on how strict you were with your gating, there may be some debris
# remaining in your preprocessed .fcs files. FlowSOM will cluster these
# together, and you may decide to exclude them from further analysis.

# Generate heatmap
plotMetaclusterMFIs(fsom_dt_dt)

# Generate UMAP
plotUMAP(fsom_dt_dt, num_cells = 5000, seed = 42)

# NOTE: In the plots panel, you may page between plots you have created using
# the arrows in the top left. This is easier than calling the functions each
# time you want to look at a previously created plot. You may also export them
# as a PDF or PNG, if you wish.

# CD8 T cells
plotLabeled2DScatter(fsom_dt, c(9,7,8,24,30,33,27,22,21,28,16,6,17,29), "metaclusters", c("BUV395-A", "Alexa Fluor 700-A"))

# CD4 T cells
plotLabeled2DScatter(fsom_dt, c(1,2,3,14,18,12,11,4,10,15,19,32,26,20), "metaclusters", c("BUV395-A", "BV650-A"))

# gdT cells
plotLabeled2DScatter(fsom_dt, c(13,45,35), "metaclusters", c("APC-Cy7-A", "BUV563-A"))

# Monocytes
plotLabeled2DScatter(fsom_dt, c(60,59,58,63,65,62,57,70,51,61,50,49), "metaclusters",  c("PE-A", "BV750-P-A"))

# B cells
plotLabeled2DScatter(fsom_dt, c(71,64,67,52,72,73,74,75), "metaclusters", c("BUV661-A", "APC-Cy7-A"))

# NK cells
plotLabeled2DScatter(fsom_dt, c(12,24,30,45,19,35,38,39,34,69,56,47,42,41), "metaclusters", c("BB700-P-A", "BUV395-A"))

# NK T cells
plotLabeled2DScatter(fsom_dt, c(12,24,30,38,39,34), "metaclusters", c("BB700-P-A", "BUV395-A"))

# cDC
plotLabeled2DScatter(fsom_dt, c(21, 31), "metaclusters", c("BV605-A", "BV750-P-A"))
plotLabeled2DScatter(fsom_dt, c(21, 31), "metaclusters", c("BV605-A", "APC-A"))

# pDC
plotLabeled2DScatter(fsom_dt, c(6, 16), "metaclusters", c("BUV737-A", "BV750-P-A"))
plotLabeled2DScatter(fsom_dt, c(6, 16), "metaclusters", c("BUV737-A", "BV605-A"))


# The editMetaclusters() function may be used to rename and merge metaclusters.
# Both of these operations are done by editing the "newLabels" parameter. For
# example, if you wanted to merge metacluster 2 into metacluster 1, and
# metacluster 9 into metacluster 4, you would call the function as follows:
#
# new_fsom_dt = UpdateMetaclusters(fsom_dt = fsom_dt, newLabels = c("1" = "2",
#                                                          "4" = "9"))
#
# Similarly, if you would like to rename metaclusters, you would do so as:
#
# new_fsom_dt = UpdateMetaclusters(fsom_dt = fsom_dt, newLabels = c("1" = "CD4 T cells",
#                                                          "4" = "CD8 T cells"))
#
# For the sake of simplicity, it is recommended that you wait until you have
# finished merging all metaclusters before you rename them, but you may wish to
# merge clusters one or a few at a time and check the results with new dimension
# reduction plots. Regardless of the approach you use, please take care to
# ensure that there are no typos in your metacluster names. If you do something
# like renaming metacluster 1 to "Monocytes" and metacluster 2 to " Monocytes"
# (note the space at the beginning of the second string), then the clusters
# will be renamed, but not merged.
#
# Note that this function may also be used to reassign individual clusters. This
# may be done by editing the "clusterAssignment" parameter, like so:
#
# new_fsom_dt = UpdateMetaclusters(fsom_dt = fsom_dt, clusterAssignment = c("54" = "NK cells",
#                                                                  "82" = "NK T cells",
#                                                                  "84" = "NK T cells))

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

# Check new plots.
plotUMAP(fsom_dt, 2000)
plotMetaclusterMFIs(fsom_dt, markers_to_cluster)

# Update names.
fsom_dt <- editTableMetaclusters(fsom_dt,
                                 new_labels = c("6" = "Monocytes",
                                                "5" = "cDC",
                                                "4" = "B cells",
                                                "17" = "NK T cells",
                                                "19" = "NK cells",
                                                "1" = "CD8 T cells",
                                                "11" = "gdT cells",
                                                "15" = "CD4 T cells",
                                                "18" = "pDC",
                                                "3" = "Undefined"),
                                 cluster_assignments = c("6" = "CD4 T cells",
                                                         "87" = "CD4 T cells",
                                                         "72" = "CD4 T cells",
                                                         "101" = "Undefined",
                                                         "7" = "Undefined",
                                                         "22" = "Undefined"))

# Final UMAP.
plotUMAP(fsom_dt, 2000)

# Get final heatmap.
plotMetaclusterMFIs(fsom_dt, markers_to_cluster)

# Generate annotated heatmap
annotateMFIHeatmap(fsom_dt, cols_to_cluster)

# Save the new table to the working directory.
# saveRDS(new_fsom_dt, paste0(dir_rds_edited, "fsom_dt_edited_26_meta_merged.rds"))
write.csv(fsom_dt, "clustered_table.csv", row.names = FALSE)

# Save the clustered .fcs files to the directory "Data/Clustered Raw/".
# SaveClustersToFCS(fsom_dt = new_fsom_dt,
#                   originalFiles = paste0(dir_prepr, prepr_files),
#                   outputDir = dir_clustr)

#########################
# DIFFERENTIAL ANALYSIS #
#########################

# Up until this section, keep clustered object as FlowSOM object?
# ^ issue: FlowSOM object doesn't track original metaclusters
# function to add item to FlowSOM object?

# Define the groups you would like to test for differential abundance and expression
comparisons <- list(
  male_vs_female = list(Sex = list("male", "female")),
  male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female")),
  ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

# Read in .csv file specifying sample information (e.g. filename, sex, age)
sample_file <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv")

# Make file compatible with flowFun functions
sample_info <- prepareSampleInfo(filepath = "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv",
                                 name.col = "Sample.Name",
                                 filename_col = "File.Name",
                                 comparisons = comparisons,
                                 samples_to_remove = NULL)

# Make design and contrasts matrices for testing
design <- makeDesignMatrix(sample_info)
contrasts <- makeContrastsMatrix(sample_info, comparisons)

# Define metaclusters to test - here, "Undefined" is excluded !!! adjust
meta_of_interest <- c("Monocytes", "CD4 T cells", "CD8 T cells", "gdT cells",
                      "B cells", "pDC", "cDC", "NK cells", "NK T cells")

# Generate matrix of sample/metacluster cell counts
counts <- makeCountMatrix(fsom_dt,
                          meta_names = meta_of_interest,
                          min_cells = 3,
                          min_samples = 4)

# Perform differential abundance analysis
da <- doDAAnalysis(design = design,
                   counts = counts,
                   contrasts = contrasts,
                   sample_df = sample_info,
                   norm_method = "TMM")

# Perform differential expression analysis
de <- doDEAnalysisNew(fsom_dt,
                      sample_df = sample_info,
                      design = design,
                      contrasts = contrasts,
                      cols_to_use = c("PHA-L", "IL10R"),
                      meta_names = meta_of_interest,
                      ctrl_input = NULL,
                      subsetted_meta = NULL,
                      save_csv = FALSE,
                      dir_tables = NULL)

################
# RECLUSTERING #
################

# If you are interested in examining a specific population in greater detail,
# you can recluster by gating on the population of interest and repeating the
# workflow outlined above. When reclustering, there is the issue of which
# markers to use for the clustering. For instance, if you decide to recluster
# T cells, markers like IgD may no longer be relevant, and not worth including.
# This script gives a few options to address this. The simplest approach is to
# use every marker in the clustering, regardless of whether or not it is relevant.
# The user may also decide to rely on prior knowledge to select which markers
# are relevant for clustering the chosen cell population. The final method is
# to perform principal component analysis (PCA) on the data, and use
# the principal components that explain the majority of the variance for
# clustering. This allows for high reproducibility while still discarding
# irrelevant information. The functions defined below implement these
# approaches.



#####################
# FUNCTIONS FOR PCA #
#####################

# The following functions may be used if you would like to use PCA for
# reclustering.


#########################
# RECLUSTERING WITH PCA #
#########################

# CD4 T CELLS

# Create a new aggregate file made up of only CD4 T cells.
agg_t = filter_frames(fsom_dt = fsom_dt,
                      dir_prepr = dir_prepr,
                      dir_clustr = dir_clustr,
                      subset_by = "metaclusters",
                      output_dir_name = "CD4 T cells B",
                      clusters_to_use = "CD4 T cells",
                      num_cells = 1000000)

# Alternatively, read in a previously created aggregate file.
agg_t = read.FCS(paste0(dir_agg, "aggregate_CD8 T cells B.fcs"))

# Perform PCA on our aggregate data.
pca = agg_pca(agg_t, markers_to_cluster)

# Plot the amount of variance explained by each principal component. In this case,
# there is not an obvious elbow, so we should look for the number of principal
# components that explain about 80% of the variance.
plot_fsom_dt_pca(pca)

# Create a new FlowSOM object using our aggregate file your chosen number
# of principal components. We will again repeat the workflow of overclustering,
# and merging clusters by examining heatmaps and dimension reduction plots.
t_fsom_dt_pca = cluster_subset_pca(aggregate = agg_t,
                                fsom_dt_filename = "fsom_dt_cd4_b",
                                prcomp = pca,
                                num_components = 11,
                                num_clus = 20,
                                xdim = 10,
                                ydim = 10,
                                dir_prepr = dir_prepr)

# Alternatively, read in a previously created FlowSOM object.
# t_fsom_dt_pca = readRDS(paste0(dir_rds_unedited, "fsom_dt_unedit_t_cells_week19.rds"))

# Plot the MST that resulted from our clustering.
t_star = PlotStars(t_fsom_dt_pca,
                   backgroundValues = t_fsom_dt_pca$metaclustering,
                   markers = markers_to_cluster)

# Plot UMAP.
t_umap = plot_umap(t_fsom_dt_pca, 2000)

# Create misc. plots of interest.
t_num = PlotNumbers(t_fsom_dt_pca, level = "metaclusters")

# Create heatmap.
t_mfi = plot_mfi_heatmap(fsom_dt = t_fsom_dt_pca,
                         markers_of_interest = markers_to_cluster)


plot_mfi_heatmap(fsom_dt = t_fsom_dt_pca,
                 markers_of_interest = c("CD45RA", "CCR7", "Fas"))

# Fas-
plotLabeled2DScatter(t_fsom_dt_pca, c(13,16), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))

# CCR7hiFas+
plotLabeled2DScatter(t_fsom_dt_pca, c(15,11), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))

# CCR7loFas+
plotLabeled2DScatter(t_fsom_dt_pca, c(12), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))

# First merging.
new_pca_fsom_dt = UpdateMetaclusters(t_fsom_dt_pca,
                                  newLabels = c("16" = "13",
                                                "15" = "11",
                                                "5" = "1",
                                                "8" = "1",
                                                "14" = "2",
                                                "20" = "19",
                                                "4" = "3",
                                                "7" = "3",
                                                "6" = "3",
                                                "10" = "3",
                                                "17" = "3"))

plot_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"))
plot_mfi_heatmap(new_pca_fsom_dt, markers_to_cluster)

# Fas-
plotLabeled2DScatter(new_pca_fsom_dt, c(13), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))
plotLabeled2DScatter(new_pca_fsom_dt, c(13), "metaclusters", c("BV750-P-A", "BUV395-A"))

# CCR7hiFas+
plotLabeled2DScatter(new_pca_fsom_dt, c(11), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))
plotLabeled2DScatter(new_pca_fsom_dt, c(11), "metaclusters", c("BV750-P-A", "BUV395-A"))

# CCR7loFas+
plotLabeled2DScatter(new_pca_fsom_dt, c(12), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))
plotLabeled2DScatter(new_pca_fsom_dt, c(12), "metaclusters", c("BV750-P-A", "BUV395-A"))

# -,-
plotLabeled2DScatter(new_pca_fsom_dt, c(3), "metaclusters", c("BV750-P-A", "BUV395-A"))

# +,-
plotLabeled2DScatter(new_pca_fsom_dt, c(19), "metaclusters", c("BV750-P-A", "BUV395-A"))

# -,+
plotLabeled2DScatter(new_pca_fsom_dt, c(1,9,18), "metaclusters", c("BV750-P-A", "BUV395-A"))

new_pca_fsom_dt = UpdateMetaclusters(new_pca_fsom_dt,
                                  newLabels = c("9" = "1",
                                                "18" = "1"),
                                  clusterAssignment = c("75" = "19",
                                                        "76" = "19",
                                                        "8" = "19",
                                                        "27" = "3",
                                                        "49" = "3",
                                                        "16" = "3",
                                                        "26" = "3",
                                                        "6" = "3"))

new_pca_fsom_dt = UpdateMetaclusters(new_pca_fsom_dt,
                                  newLabels = c("11" = "CD45RA+, CCR7+ (CCR7hiFas+)",
                                                "12" = "CD45RA+, CCR7+ (CCR7loFas+)",
                                                "13" = "CD45RA+, CCR7+ (Fas-)",
                                                "1" = "CD45RA-, CCR7+",
                                                "2" = "CD45RA-, CCR7+",
                                                "19" = "CD45RA+, CCR7-",
                                                "3" = "CD45RA-, CCR7-"))

saveRDS(new_pca_fsom_dt, paste0(dir_rds_edited, "fsom_dt_cd4_edited_B.rds"))

# Create final heatmap.
anno_hm = add_merging_anno(t_mfi, t_fsom_dt_pca, new_pca_fsom_dt)

# Save final clustering to "RDS/Edited/" directory.
saveRDS(new_pca_fsom_dt, paste0(dir_rds_edited, "fsom_dt_t_edited_week19.rds"))

##############
# CD8 T cells

agg_cd8 = filter_frames(fsom_dt = fsom_dt,
                        dir_prepr = dir_prepr,
                        dir_clustr = dir_clustr,
                        subset_by = "metaclusters",
                        output_dir_name = "CD8 T cells B",
                        clusters_to_use = "CD8 T cells",
                        num_cells = 1000000)

# Alternatively, read in a previously created aggregate file.
# agg_t = read.FCS(paste0(dir_agg, "aggregate_T_cells_week19.fcs"))

# Perform PCA on our aggregate data.
pca = agg_pca(agg_cd8, markers_to_cluster)

# Plot the amount of variance explained by each principal component. In this case,
# there is not an obvious elbow, so we should look for the number of principal
# components that explain about 80% of the variance.
plot_fsom_dt_pca(pca)

# Create a new FlowSOM object using our aggregate file your chosen number
# of principal components. We will again repeat the workflow of overclustering,
# and merging clusters by examining heatmaps and dimension reduction plots.
cd8_fsom_dt_pca = cluster_subset_pca(aggregate = agg_cd8,
                                  fsom_dt_filename = "fsom_dt_cd8_b",
                                  prcomp = pca,
                                  num_components = 10,
                                  num_clus = 20,
                                  xdim = 10,
                                  ydim = 10,
                                  dir_prepr = dir_prepr)

# Plot the MST that resulted from our clustering.
cd8_star = PlotStars(cd8_fsom_dt_pca,
                     backgroundValues = cd8_fsom_dt_pca$metaclustering,
                     markers = markers_to_cluster)

# Plot UMAP.
cd8_umap = plot_umap(cd8_fsom_dt_pca, 2000)

# Create heatmap.
cd8_mfi = plot_mfi_heatmap(fsom_dt = cd8_fsom_dt_pca,
                           markers_of_interest = markers_to_cluster)

plot_mfi_heatmap(cd8_fsom_dt_pca, c("CD45RA", "CCR7", "Fas"))

# First merging.
new_pca_fsom_dt = UpdateMetaclusters(cd8_fsom_dt_pca,
                                  newLabels = c("6" = "3",
                                                "4" = "3",
                                                "16" = "11",
                                                "15" = "11",
                                                "13" = "2",
                                                "8" = "2",
                                                "20" = "19",
                                                "14" = "12"))

plot_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"))

# Fas-
plotLabeled2DScatter(new_pca_fsom_dt, c(5), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))
plotLabeled2DScatter(new_pca_fsom_dt, c(5), "metaclusters", c("BV750-P-A", "BUV395-A"))

# CCR7hiFas+
plotLabeled2DScatter(new_pca_fsom_dt, c(3), "metaclusters", c("BV750-P-A", "PE-Cy7-A"))
plotLabeled2DScatter(new_pca_fsom_dt, c(3), "metaclusters", c("BV750-P-A", "BUV395-A"))

# CCR7loFas+
plotLabeled2DScatter(new_pca_fsom_dt, c(9), "metaclusters", c("BV750-P-A", "PE-Cy7-A"), plot_meta_labels = FALSE)
plotLabeled2DScatter(new_pca_fsom_dt, c(9), "metaclusters", c("BV750-P-A", "BUV395-A"), plot_meta_labels = FALSE)

plotLabeled2DScatterer_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"), c(3,5))

new_pca_fsom_dt = UpdateMetaclusters(new_pca_fsom_dt,
                                  newLabels = c("12" = "7",
                                                "19" = "7",
                                                "17" = "2",
                                                "11" = "9"),
                                  clusterAssignment = c("28" = "CD45RA+, CCR7+ (CCR7loFas+)",
                                                        "37" = "CD45RA+, CCR7+ (CCR7loFas+)"))

plot_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"))

plotLabeled2DScatterer_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"), c(1,7))

new_pca_fsom_dt = UpdateMetaclusters(new_pca_fsom_dt,
                                  newLabels = c("10" = "2",
                                                "18" = "2"),
                                  clusterAssignment = c("40" = "3",
                                                        "23" = "2",
                                                        "47" = "2"))

plot_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"))

plotLabeled2DScatter(new_pca_fsom_dt, c(1,7), "metaclusters", c("BV750-P-A", "PE-Cy7-A"), plot_meta_labels = FALSE)
plotLabeled2DScatter(new_pca_fsom_dt, c(1,7), "metaclusters", c("BV750-P-A", "BUV395-A"), plot_meta_labels = FALSE)

new_pca_fsom_dt = UpdateMetaclusters(new_pca_fsom_dt,
                                  newLabels = c("3" = "CD45RA+, CCR7+ (CCR7hiFas+)",
                                                "9" = "CD45RA-, CCR7+",
                                                "5" = "CD45RA+, CCR7+, Fas-",
                                                "2" = "CD45RA+, CCR7-",
                                                "1" = "CD45RA-, CCR7-",
                                                "7" = "CD45RA-, CCR7-"))

plot_mfi_heatmap(new_pca_fsom_dt, c("CD45RA", "CCR7", "Fas"))

anno_cd8_hm = add_merging_anno(cd8_mfi, cd8_fsom_dt_pca, new_pca_fsom_dt)

# Save final clustering to "RDS/Edited/" directory.
saveRDS(new_pca_fsom_dt, paste0(dir_rds_edited, "fsom_dt_cd8_edited_B.rds"))

############################
# RECLUSTERING WITHOUT PCA #
############################

# THIS IS JUST AN EXAMPLE, NOT ACTUAL ANALYSIS

# The process for backgating without PCA is practically identical to what we
# have above, except for there is no need to create a prcomp object, and we use
# cluster_subset() rather than cluster_subset_pca(). We may use the same
# aggregate file we created above with the filter_frames() function.
b_fsom_dt = cluster_subset(aggregate = agg_b,
                        fsom_dt_filename = "fsom_dt_b_no_pca_unedited",
                        markers_to_cluster= markers_to_cluster,
                        num_clus = 10)

# Alternatively, read in a previously created FlowSOM object.
# cd8_fsom_dt = readRDS(paste0(dir_rds_unedited, "fsom_dt_cd8_20clus_unedited.rds"))

# Plot resulting MST.
PlotStars(b_fsom_dt, backgroundValues = b_fsom_dt$metaclustering)

# Create plots to reference for merging.
plot_umap(b_fsom_dt, 2000)
plot_mfi_heatmap(b_fsom_dt, markers_to_cluster)

# First merging.
new_cd8_fsom_dt = UpdateMetaclusters(cd8_fsom_dt,
                                  newLabels = c("13" = "10",
                                                "2" = "1",
                                                "19" = "15",
                                                "16" = "15"))


plot_mfi_heatmap(new_cd8_fsom_dt, c(markers_to_cluster, "PHA-L", "Fas", "IL10R"))
plot_mfi_heatmap(new_cd8_fsom_dt, c("CD45RA", "CCR7", "Fas"))
plot_umap(new_cd8_fsom_dt, 2000)

# Second merging.
new_cd8_fsom_dt = UpdateMetaclusters(new_cd8_fsom_dt,
                                  newLabels = c("14" = "11",
                                                "9" = "4",
                                                "7" = "1",
                                                "6" = "3",
                                                "8" = "5"))

plot_mfi_heatmap(new_cd8_fsom_dt, c("CD45RA", "CCR7", "Fas"))
plot_umap(new_cd8_fsom_dt, 2000)

# Third merging.
new_cd8_fsom_dt = UpdateMetaclusters(new_cd8_fsom_dt,
                                  newLabels = c("3" = "1",
                                                "11" = "4"))

plot_mfi_heatmap(new_cd8_fsom_dt, c("CD45RA", "CCR7", "Fas"))
plot_umap(new_cd8_fsom_dt, 2000)

# Plot MST of final clustering.
PlotStars(new_cd8_fsom_dt, backgroundValues = new_cd8_fsom_dt$metaclustering)

# Save final clustering to "RDS/Edited/" directory.
saveRDS(new_pca_fsom_dt, paste0(dir_rds_edited, "fsom_dt_cd8_edited.rds"))


###################
# MISC. FUNCTIONS #
###################

# plot_prcntgs() produces a bar plot that displays cell type proportions by
# sample. To use it, you must pass your FlowSOM object, and the name of the
# directory containing your preprocessed files. By default, the names used for
# each sample will be its file name. However, you may rename the samples by
# providing a regular expression with a capturing group, which will be applied
# to each file name.
plot_prcntgs = function(fsom_dt, dir_prepr, reg_expr = NULL) {
  sample_names = list.files(dir_prepr)
  sample_prcntgs = c()
  for (i in 1:length(sample_names)) {
    print(paste0("Processing ", sample_names[i], "..."))
    ind = which(fsom_dt$data[, "File"] == i)
    fsom_dt_subset = FlowSOMSubset(fsom_dt, ind)
    counts = GetCounts(fsom_dt_subset, level = "metaclusters")
    sample_prcntgs = rbind(sample_prcntgs, counts/sum(counts))
  }

  if (!is.null(reg_expr)) {
    sample_names = sub(reg_expr, "\\1", sample_names)
  }

  sample_prcntgs = cbind(sample_names, sample_prcntgs)
  sample_prcntgs = as.data.frame(sample_prcntgs)

  prcntgs_long = tidyr::pivot_longer(sample_prcntgs,
                                     cols = -sample_names,
                                     names_to = "Cell Type",
                                     values_to = "Proportion")
  prcntgs_long$Proportion = as(prcntgs_long$Proportion, "numeric")

  ggplot(prcntgs_long, aes(x = sample_names, y = Proportion, fill = `Cell Type`)) +
    geom_bar(stat = "identity") +
    labs(x = "Sample Names", y = "Proportion") +
    ggtitle("Cell Types by Sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    scale_fill_viridis_d(option = "turbo")
}

plot_prcntgs(new_fsom_dt, dir_prepr, "LN_PHA (.*?)\\.fcs")


# remove_filepath() is a helper function that removes everything but the file
# name from a file path.
remove_filepath = function(filenames) {
  contains_path = grepl("/", filenames)
  new_filenames = ifelse(contains_path, sub(".*/", "", filenames), filenames)
  return(new_filenames)
}

# plot_group_dist() plots group proportions by cluster. For each cluster, a pie
# chart displaying the proportions of each given group within said cluster is
# created. They may be visualized as either an MST or a grid. The pie chart in
# the bottom left may serve as a point of reference, as it shows what we would
# expect the proportions to be if there was no difference between groups.
#
# The parameters represent the following:
# fsom_dt: The FlowSOM object whose file distribution you would like to plot.
#
# directory: The directory containing the files that were used for clustering.
#
# groups: A list of lists, specifying the groups of interest and their file
#   names.
#
# view: A string specifying which version of the plot you would like returned _
#   either "MST" or "grid".
plot_group_dist = function(fsom_dt, directory, groups, view) {
  file_names = list.files(path = directory,
                          full.names = FALSE)
  all_group_files = c()
  for (group in groups) {
    all_group_files = c(all_group_files, group)
  }

  sample_names = file_names[fsom_dt$data[,"File"]]
  ind = which(paste0(directory, sample_names) %in% all_group_files)

  if (length(ind) < length(sample_names)) {
    sample_names = sample_names[ind]
    fsom_dt = FlowSOMSubset(fsom_dt, ind)
  }

  group_indicator = rep("Unknown", length(sample_names))

  # check names for both "file_names" and "groups"
  for (i in 1:length(groups)) {
    groups[[i]] = remove_filepath(groups[[i]])
    group_indicator[sample_names %in% unlist(groups[names(groups)[i]])] = names(groups)[i]
  }

  file_pie_plot = PlotPies(fsom_dt = fsom_dt,
                           cellTypes = factor(group_indicator),
                           colorPalette = viridis(length(groups), option = "inferno"),
                           equalNodeSize = TRUE,
                           maxNodeSize = ifelse(view == "MST", 0.7, 1),
                           view = view)

  group_proportions = sapply(groups, length)
  group_proportions = group_proportions / sum(group_proportions)

  angles = cumsum(2 * pi * group_proportions)
  start_angles = c(0, head(angles, -1))

  file_pie_plot = AddStarsPies(p = file_pie_plot,
                               arcs = data.frame(x0 = 0,
                                                 y0 = 0,
                                                 start = start_angles,
                                                 end = angles,
                                                 value = 1,
                                                 Markers = names(groups)),
                               colorPalette = viridis(length(groups), option = "inferno")
  )

  return(file_pie_plot)
}

# Plot Ctrl/MIBC groups
ctrl_files = list.files(path = dir_prepr,
                        pattern = "Bld_PHA Naive.*\\.fcs",
                        full.names = TRUE)
bbn_files = list.files(path = dir_prepr,
                       pattern = "Bld_PHA BBN.*\\.fcs",
                       full.names = TRUE)
disease_groups = list("Naive" = ctrl_files,
                      "BBN"= bbn_files)

plot_group_dist(new_fsom_dt, dir_prepr, disease_groups, "MST")

# Plot male/female groups
file_info = read.csv("Info/mouse_week19_sample_info.csv", check.names = FALSE)
female_group = paste0(dir_prepr, file_info[file_info$Sex == "Female", ]$`Filename`)
male_group = paste0(dir_prepr, file_info[file_info$Sex == "Male", ]$`Filename`)

sex_groups = list("Female" = female_group,
                  "Male"= male_group)

plot_group_dist(new_fsom_dt, dir_prepr, sex_groups, "MST")
plot_group_dist(new_pca_fsom_dt, dir_prepr, sex_groups, "MST")


# plot_file_dist() plots file distribution per cluster. Like the above function,
# the result may be visualized as either an MST or a grid. The pie chart in
# the bottom left may serve as a point of reference, as it shows what we would
# expect the proportions to be if there was no difference between samples, and
# if the clustering went well. Note that if a metacluster or cluster is made up
# mostly of one sample's cells, this does not necessarily mean that the
# clustering is flawed, and may instead tell us that the relevant sample has
# unique biological characteristics.
#
# To use the function, the parameters fsom_dt, directory, and view, work identical
# to the plot_group_dist() function defined above. However, the function may
# take in an additional parameter, "reg_expr", which may be used to specify a
# regular expression that will be used to rename the files for the legend. If no
# regular expression is provided, the full file names will be used instead.
plot_file_dist = function(fsom_dt, directory, view, reg_expr = NULL) {
  sample_names = list.files(path = directory,
                            full.names = FALSE)

  if (!is.null(reg_expr)) {
    sample_names = sub(reg_expr, "\\1", sample_names)
  }

  file_num = length(sample_names)

  file_pie_plot = PlotPies(fsom_dt = fsom_dt,
                           cellTypes = factor(sample_names[fsom_dt$data[,"File"]]),
                           colorPalette = viridis(file_num, option = "turbo"),
                           equalNodeSize = TRUE,
                           maxNodeSize = ifelse(view == "MST", 0.7, 1),
                           view = view)

  proportions = sapply(1:length(sample_names), function(i) {
    num_cells = length(which(fsom_dt$data[, "File"] == i))
    return(num_cells)
  })

  proportions = proportions / sum(proportions)

  angles = cumsum(2 * pi * proportions)
  start_angles = c(0, head(angles, -1))

  file_pie_plot = AddStarsPies(p = file_pie_plot,
                               arcs = data.frame(x0 = 0,
                                                 y0 = 0,
                                                 start = start_angles,
                                                 end = angles,
                                                 value = 1,
                                                 Markers = sample_names),
                               colorPalette = viridis(file_num, option = "turbo"))

  return(file_pie_plot)
}

plot_file_dist(new_fsom_dt, dir_prepr, "grid", "Bld_PHA (.*?)\\.fcs")
plot_file_dist(new_pca_fsom_dt, dir_prepr, "grid", "Bld_PHA (.*?)\\.fcs")


# Helper for extra plotting functions below.
get_density = function(x, y, ...) {
  dens = MASS::kde2d(x, y, ...)
  ix = findInterval(x, dens$x)
  iy = findInterval(y, dens$y)
  ii = cbind(ix, iy)
  return(dens$z[ii])
}

# Plot t-SNEs by group. The parameters represent the following:
#
# fsom_dt: The FlowSOM object whose result you would like to plot.
#
# clustered_markers: The markers that were used for clustering. In most scenarios
#   this will be your `markers_to_cluster` variable defined at the beginning of
#   the file.
#
# prepr_files: All preprocessed files that were clustered. Most often this should
#   be the `prepr_files` variable defined at the beginning of the file.
#
# groups: A list specifying the groups of interest.
#
# color_by: A string specifying how you would like the t-SNE plot to be colored.
#   Your options are "metacluster", "density", or a marker or channel of interest.
#
# num_cells: The number of cells you would like to be sampled for each group.
plot_group_tsne = function(fsom_dt, clustered_markers, prepr_files, groups,
                           color_by, num_cells) {
  set.seed(42)

  all_ind = c()

  for (group in names(groups)) {
    groups[[group]] = remove_filepath(groups[[group]])
    group_ind = match(groups[[group]], prepr_files)
    temp = which(fsom_dt$data[, "File"] %in% group_ind)
    samp = sample(temp, num_cells)
    all_ind = c(all_ind, samp)
  }

  dat = fsom_dt$data[all_ind, GetChannels(fsom_dt, clustered_markers)]

  tsne = Rtsne(dat, seed = 42)

  group_vec = rep(1, num_cells*length(groups))
  Density = c()
  for (i in 1:length(groups)) {
    ind = ((num_cells*(i-1))+1):(num_cells*i)
    group_vec[ind] = names(groups)[i]
    if (color_by == "density") {
      Density = c(Density, get_density(tsne$Y[ind,1], tsne$Y[ind,2], n = 100))
    }
  }

  if (color_by == "metacluster") {
    meta_vec = GetMetaclusters(fsom_dt)[all_ind]
    tsne_df = data.frame(tsne$Y, group = group_vec, Metacluster = meta_vec)

    ggplot(tsne_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster),
                                                    pointsize = 4) +
      facet_wrap(~group) +
      theme_minimal()
  } else if (color_by == "density") {
    tsne_df = data.frame(tsne$Y, group = group_vec)
    # Density = get_density(tsne$Y[,1], tsne$Y[,2], n = 100)

    ggplot(tsne_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Density),
                                                    pointsize = 4) +
      scale_color_viridis(option = "H") +
      facet_wrap(~group) +
      theme_minimal()
  } else {
    if (color_by %in% colnames(fsom_dt$data)) {
      channel = color_by
    } else if (color_by %in% GetMarkers(fsom_dt, colnames(fsom_dt$data))) {
      channel = GetChannels(fsom_dt, color_by)
    } else {
      stop("The value given to parameter `color_by` is invalid.")
    }
    marker_vec = fsom_dt$data[all_ind, channel]
    tsne_df = data.frame(tsne$Y, group = group_vec, Expression = marker_vec)

    ggplot(tsne_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Expression),
                                                    pointsize = 4) +
      scale_color_viridis(option = "H") +
      facet_wrap(~group) +
      theme_minimal()
  }
}

# Plot UMAPs by group. The parameters represent the following:
#
# fsom_dt: The FlowSOM object whose result you would like to plot.
#
# clustered_markers: The markers that were used for clustering. In most scenarios
#   this will be your `markers_to_cluster` variable defined at the beginning of
#   the file.
#
# prepr_files: All preprocessed files that were clustered. Most often this should
#   be the `prepr_files` variable defined at the beginning of the file.
#
# groups: A list specifying the groups of interest, in terms of file names.
#
# color_by: A string specifying how you would like the t-SNE plot to be colored.
#   Your options are "metacluster", "density", or a marker or channel of interest.
#
# num_cells: The number of cells you would like to be sampled for each group.
plot_group_umap = function(fsom_dt, clustered_markers, prepr_files, groups,
                           color_by, num_cells, layout = NULL) {
  set.seed(42)

  all_ind = c()

  for (group in names(groups)) {
    groups[[group]] = remove_filepath(groups[[group]])
    group_ind = match(groups[[group]], prepr_files)
    temp = which(fsom_dt$data[, "File"] %in% group_ind)
    samp = sample(temp, num_cells)
    all_ind = c(all_ind, samp)
  }

  dat = fsom_dt$data[all_ind, GetChannels(fsom_dt, clustered_markers)]

  if (!is.null(layout)) {
    umap = umap(dat, init = layout)
    print("A")
  } else {
    umap = umap(dat)
  }


  group_vec = rep(1, num_cells*length(groups))
  Density = c()
  for (i in 1:length(groups)) {
    ind = ((num_cells*(i-1))+1):(num_cells*i)
    group_vec[ind] = names(groups)[i]
    if (color_by == "density") {
      Density = c(Density, get_density(umap$layout[ind,1], umap$layout[ind,2], n = 100))
    }
  }

  if (color_by == "metacluster") {
    meta_vec = GetMetaclusters(fsom_dt)[all_ind]
    umap_df = data.frame(umap$layout, group = group_vec, Metacluster = meta_vec)

    print(ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster),
                                                          pointsize = 4) +
            facet_wrap(~group) +
            theme_minimal())
  } else if (color_by == "density") {
    umap_df = data.frame(umap$layout, group = group_vec)
    #Density = get_density(umap$layout[,1], umap$layout[,2], n = 100)

    ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Density),
                                                    pointsize = 4) +
      scale_color_viridis(option = "H") +
      facet_wrap(~group) +
      theme_minimal()
  } else {
    if (color_by %in% colnames(fsom_dt$data)) {
      channel = color_by
    } else if (color_by %in% GetMarkers(fsom_dt, colnames(fsom_dt$data))) {
      channel = GetChannels(fsom_dt, color_by)
    } else {
      stop("The value given to parameter `color_by` is invalid.")
    }
    marker_vec = fsom_dt$data[all_ind, channel]
    umap_df = data.frame(umap$layout, group = group_vec, Expression = marker_vec)

    ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Expression),
                                                    pointsize = 4) +
      scale_color_viridis(option = "H") +
      facet_wrap(~group) +
      theme_minimal()
  }
  return(umap)
}

plot_group_umap(fsom_dt_cd8, markers_to_cluster, prepr_files, sex_groups, "PHA-L", 2000)

# Automatically gate on a population of interest. Must specify two channels
# and whether the population is positive or negative for each.
subset_by_marker = function(fsom_dt, aggregate, metaclusters, channels, position) {
  inds = which(GetMetaclusters(fsom_dt) %in% metaclusters)
  aggregate = aggregate[inds, ]
  pop = flowDensity(aggregate, channels, position)
  new = aggregate[pop@index, ]
  print(paste0(round(pop@proportion, 2), "% of cells in the given metaclusters",
               " are in the population of interest."))
  return(pop)
}


test = subset_by_marker(fsom_dt_cd8, agg_cd8, c("CD45RA-, CCR7-"), c("Alexa Fluor 700-A", "PE-A"),
                        c(TRUE, TRUE))
p = plotDens(test, c("Alexa Fluor 700-A", "PE-A"))


# Plot gates resulting from above function.
plot_flowdens_gates = function(fsom_dt, aggregate, metaclusters, flowdens_obj, channel1, channel2, ncells) {
  aggregate = aggregate[which(GetMetaclusters(fsom_dt) %in% metaclusters), ]
  df = data.frame(x = exprs(aggregate)[, channel1],
                  y = exprs(aggregate)[, channel2])
  if (ncells > nrow(df)) {
    i = 1:nrow(df)
  } else {
    i = sample(nrow(df), ncells)
  }

  df = df[i, ]

  dens = get_density(df$x, df$y, n = 100)

  plot = ggplot(df, aes(x = x, y = y, color = dens)) +
    geom_point(size = 0.5) +
    xlab(get_marker(aggregate, channel1)) +
    ylab(get_marker(aggregate, channel2)) +
    geom_hline(mapping = aes(yintercept = flowdens_obj@gates[1])) +
    geom_vline(mapping = aes(xintercept = flowdens_obj@gates[2])) +
    scale_color_viridis(option = "H") +
    theme_minimal() + theme(legend.position = "none")
  ggpubr::ggarrange(plot)
}

plot_flowdens_gates(fsom_dt_cd8, agg_cd8, c("CD45RA-, CCR7-"), test, "Alexa Fluor 700-A", "PE-A", 5000)

# Alternatively, you may draw a gate manually. Click two diagonal points to draw
# a rectangular gate.
subset_by_marker_manual = function(fsom_dt, aggregate, metaclusters, channels) {
  inds = which(GetMetaclusters(fsom_dt) %in% metaclusters)
  aggregate = aggregate[inds, ]

  polygon_gate = CytoExploreR::cyto_gate_draw(x = aggregate,
                                              alias = "population",
                                              type = "rectangle",
                                              channels = channels)
  polygon_gate = polygon_gate$`population`

  new = aggregate[flowCore::filter(aggregate, polygon_gate)@subSet, ]

  print(paste0(round(nrow(new)/nrow(aggregate)*100, 2), "% of cells in the given metacluster(s)",
               " are in the population of interest."))

  return(new)
}

test2 = subset_by_marker_manual(fsom_dt_cd8, agg_cd8, c("CD45RA-, CCR7-"), c("Alexa Fluor 700-A", "PE-A"))

#####################
plot_group_umap = function(fsom_dt, clustered_markers, factors, num_cells) {
  set.seed(42)

  groups = levels(factors$group)

  all_ind = c()

  for (group in groups) {
    grp_samples = which(factors[, "group"] == group)
    inds = which(fsom_dt$data[, "File"] %in% grp_samples)
    samp = sample(inds, num_cells)
    all_ind = c(all_ind, samp)
  }

  dat = fsom_dt$data[all_ind, GetChannels(fsom_dt, clustered_markers)]

  umap = umap(dat)

  meta_vec = GetMetaclusters(fsom_dt)[all_ind]
  meta_vec = factor(meta_vec, levels = c("CD4 T cells", "CD8 T cells", "gdT cells",
                                         "B cells", "NK cells", "NK T cells", "Monocytes",
                                         "cDC", "pDC", "Undefined"))
  umap_df = data.frame(umap$layout, Metacluster = meta_vec, Indices = all_ind)

  print(ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster),
                                                        pointsize = 2) +
          theme_void())

  return(umap_df)
}

getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#####
plotGroupUMAPs <- function(fsom_dt, sample_df, grps_of_interest, umap = NULL,
                           color_by = "density", num_cells = NULL, seed = NULL) {
  # Set seed if desired
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # If no UMAP given
  if (is.null(umap)) {
    fsom_dt_inds <- c()
    group_vec <- rep(1, length(grps_of_interest) * num_cells)
    for (group in names(grps_of_interest)) {
      # Get sample indices belonging to current group
      # Note that `rownames(sample_df)` corresponds to the sample's number in `fsom_dt$data[, "File"]`
      grp_samples <- rownames(sample_df)[which(sample_df$group %in% grps_of_interest[[group]])]
      # Get cell indices of `fsom_dt` belonging to current group
      inds <- which(fsom_dt$data[, "File"] %in% grp_samples)
      # Sample `num_cells` for current group of interest
      samp <- sample(inds, num_cells)
      fsom_dt_inds <- c(fsom_dt_inds, samp)

      # Assign group name for each cell
      upper_bound <- which(names(grps_of_interest) == group) * num_cells
      grp_inds <- seq((upper_bound - num_cells + 1), upper_bound)
      group_vec[grp_inds] <- group
    }

    # Subset data and generate parent UMAP
    fsom_dt <- FlowSOM::FlowSOMSubset(fsom_dt, fsom_dt_inds)
    umap <- plotUMAP(fsom_dt, length(fsom_dt_inds))
    umap_df <- umap$data

  } else { # if a parent UMAP was given
    umap_df <- umap$data

    # Get group name for each cell
    group_vec <- fsom_dt$data[umap_df$Indices, "File"]
    for (group in names(grps_of_interest)) {
      grp_samples <- rownames(sample_df)[which(sample_df$group %in% grps_of_interest[[group]])] # samples belonging to current group
      inds <- which(group_vec %in% grp_samples) # cell indices of `dat` belonging to current group
      group_vec[inds] <- group

    }

    # Resample cells if necessary
    tab <- table(factor(group_vec, levels = names(grps_of_interest)))
    grp_inds <- c()

    # Reduce number of cells such that they are equal across groups
    # Reduce number of cells such that they are equal across groups
    if (any(tab == 0)) {
      stop(paste("There are no cells in group", names(tab)[which(tab == 0)],
                 "in the given UMAP. Generate a new one by either setting `umap = NULL`,",
                 "or using the function `plotUMAP()`."))
    } else if (all(tab > num_cells)) { # if all groups have more cells than `num_cells`
      subtra <- num_cells
      grp_inds <- 1:length(grps_of_interest)
    } else if (any(tab > min(tab))) { # if groups have disproportionate cell counts
      subtra <- min(tab)
      grp_inds <- which(tab > min(tab))

      print(paste("There are some groups with less than", num_cells,
                  "cells. UMAPs were plotted with", subtra, "cells each instead."))
    }
    for (grp_ind in grp_inds) {
      diff <- tab[grp_ind] - subtra
      old_inds <- which(group_vec == names(tab[grp_ind]))
      samp <- sample(old_inds, diff) # indices of points to remove

      umap_df <- umap_df[-samp, ]
      group_vec <- group_vec[-samp]
    }
  }

  # Initialize parameters
  X1 <- X2 <- values <- NULL

  # Get point values and legend name
  if (color_by == "density") {
    # Calculate density values
    Density <- rep(1, length(group_vec))
    for (group in names(grps_of_interest)) {
      inds <- which(group_vec == group)
      Density[inds] <- getDensity(umap_df[inds, 1], umap_df[inds, 2], n = 100)
    }

    # Append density column to data frame
    umap_df <- data.frame(umap_df, group = group_vec, values = Density/max(Density))
    legend_name <- "Density"

  } else {
    if (color_by %in% colnames(fsom_dt$data)) {
      channel <- color_by
    } else if (color_by %in% FlowSOM::GetMarkers(fsom_dt, colnames(fsom_dt$data))) {
      channel <- FlowSOM::GetChannels(fsom_dt, color_by)
    } else {
      stop("The value given to parameter `color_by` is invalid.")
    }

    # Get expression data and append column to data frame
    marker_vec <- fsom_dt$data[umap_df$Indices, channel]
    umap_df <- data.frame(umap_df, group = group_vec, values = marker_vec)
    legend_name <- channel

  }

  # Draw plot
  p <- ggplot2::ggplot(umap_df) +
    scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = values),
                                  pointsize = 2) +
    ggplot2::facet_wrap(~group) +
    viridis::scale_color_viridis(option = "H", name = legend_name) +
    ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio = 1, strip.text = ggplot2::element_text(size = 12))

  return(p)
}

# For full clustering with 26 metaclusters.
full_umap_b = plot_group_umap(new_fsom_dt, markers_to_cluster, factors, 7500)
init_full_umap = plot_group_umap(new_fsom_dt, markers_to_cluster, factors, 7500)

p = ggplot(full_umap_b) +
  #scale_color_viridis(discrete = TRUE, option = "H", alpha = 0.5) +
  scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster),
                                pointsize = 2) +
  #scale_fill_viridis(discrete = TRUE, option = "G", begin = 0.2) +
  theme_void() + theme(aspect.ratio = 1)

# For CD8 T cell reclustering.
full_umap = plot_group_umap(new_fsom_dt, markers_to_cluster, factors, 7500)

plot_group_umaps = function(fsom_dt, umap_df, factors, grps_of_interest, channel, num_cells) {
  groups = levels(factors$group)
  cells_per_group = nrow(umap_df)/length(groups)

  all_inds = c()

  for (vec in grps_of_interest) {
    all_grp_inds = c()
    for (group in vec) {
      grp_ind = which(groups == group)
      upper_bound = grp_ind * cells_per_group

      all_grp_inds = c(all_grp_inds, seq((upper_bound-cells_per_group+1), upper_bound))
    }
    all_inds = c(all_inds, sample(all_grp_inds, num_cells))
  }

  umap_df = umap_df[all_inds, ]

  group_vec = rep(1, num_cells*length(grps_of_interest))
  for (i in 1:length(grps_of_interest)) {
    ind = ((num_cells*(i-1))+1):(num_cells*i)
    group_vec[ind] = names(grps_of_interest)[i]
  }

  marker_vec = fsom_dt$data[umap_df$Indices, channel]
  umap_df = data.frame(umap_df, group = group_vec, Expression = marker_vec)

  ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Expression),
                                                  pointsize = 2) +
    scale_color_viridis(option = "H") +
    facet_wrap(~group) +
    theme_void() +
    theme(aspect.ratio = 1, strip.text = element_text(size = 12))
}

grps_of_interest = list("Female" = c("female_Ctrl_X", "female_MIBC_No.NAC", "female_MIBC_NAC"),
                        "Male" = c("male_Ctrl_X", "male_MIBC_No.NAC", "male_MIBC_NAC"))

mibc_sex_grps_of_interest = list("Female MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC"),
                                 "Male MIBC" = c("male_MIBC_No.NAC", "male_MIBC_NAC"))

sex_four_group = list("Male Ctrl" = c("male_Ctrl_X"),
                      "Female Ctrl" = c("female_Ctrl_X"),
                      "Male MIBC" = c("male_MIBC_No.NAC", "male_MIBC_NAC"),
                      "Female MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC"))

dis_grps_of_interest = list("Control" = c("female_Ctrl_X", "male_Ctrl_X"),
                            "MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC", "male_MIBC_No.NAC", "male_MIBC_NAC"))

plot_group_umaps(new_fsom_dt, full_umap_b, factors, grps_of_interest, "FITC-A", 22500)
plot_group_umaps(new_fsom_dt, full_umap_b, factors, dis_grps_of_interest, "FITC-A", 15000)
plot_group_umaps(new_fsom_dt, full_umap_b, factors, mibc_sex_grps_of_interest, "FITC-A", 15000)

pdf("Analysis Results/Misc Plots/UMAPs/full_plots_merged.pdf")
# Colored by metacluster
ggplot(init_full_umap) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster),
                                                       pointsize = 2) +
  theme_void() + theme(aspect.ratio = 1, strip.text = element_text(size = 12))

# Grouped by disease.
plot_group_umaps(new_fsom_dt, full_umap_b, factors, dis_grps_of_interest, "FITC-A", 15000)

# Grouped by sex.
plot_group_umaps(new_fsom_dt, full_umap_b, factors, grps_of_interest, "FITC-A", 22500)

# Grouped by sex with MIBC.
plot_group_umaps(new_fsom_dt, full_umap_b, factors, mibc_sex_grps_of_interest, "FITC-A", 15000) + theme(strip.text = element_text(size = 12))

# Bar plot, grid lines.
bar_plot_fun("Analysis Results/cd8_pha_dMFI_linear_table.csv", factors, sex_four_group) +
  xlab("Cell Type") + ylab("dMFI") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 0.8, size = 8), aspect.ratio = 0.6) +
  scale_shape_manual(values = c(16,17,15,18))

# Bar plot, no grid lines.
bar_plot_fun("Analysis Results/cd8_pha_dMFI_linear_table.csv", factors, sex_four_group) +
  xlab("Cell Type") + ylab("dMFI") + theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8), aspect.ratio = 0.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8500)) +
  scale_shape_manual(values = c(16,17,15,18))
dev.off()

##################
# bar plots

se_fun = function(x) {
  sd = sd(x)
  se = sd/sqrt(length(x))
  return(se)
}

bar_plot_fun = function(csv_file, factors, grps_of_interest) {
  orig_df = read.csv(csv_file, check.names = FALSE)
  orig_df = cbind(orig_df, factor_group = factors$group)

  df = orig_df[,-1]

  for (col in colnames(df)) {
    if (col != "factor_group") {
      df[, col] = as.numeric(df[, col])
    }
  }

  # Ensure "group" column is treated as a factor.
  df = df %>%
    mutate(factor_group = as.factor(factor_group))

  df$Group = factor(seq(1:nrow(df)), levels = names(grps_of_interest))
  for (i in 1:nrow(df)) {
    group = as.character(df$factor_group[i])
    for (j in 1:length(grps_of_interest)) {
      if (group %in% grps_of_interest[[j]]) {
        #print(names(grps_of_interest)[j])
        df$Group[i] = names(grps_of_interest)[j]
      }
    }
  }
  #print(head(df))

  df = df[, -(which(colnames(df) == "factor_group"))]
  point_df = pivot_longer(df, cols = -Group)

  # new_df = df %>%
  #   group_by(Group) %>%
  #   summarise_if(is.numeric, list(Median = median, SD = sd), na.rm = TRUE) %>%
  #   ungroup()
  median_df = df %>%
    group_by(Group) %>%
    summarise_if(is.numeric, median, na.rm = TRUE) %>%
    ungroup()

  se_df = df %>%
    group_by(Group) %>%
    summarise_if(is.numeric, se_fun) %>%
    ungroup()

  med_dat_long = pivot_longer(median_df, cols = -Group)
  se_dat_long = pivot_longer(se_df, cols = -Group)

  upper_lim = max(med_dat_long$value)
  upper_lim = upper_lim + (upper_lim/2)

  ggplot(med_dat_long, aes(x = name, y = value, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_point(data = point_df,
               mapping = aes(name, value, shape = Group, fill = Group),
               position = position_jitterdodge(jitter.width = 0.05, dodge.width = 1, seed = 42),
               size = 1,
               inherit.aes = FALSE) +
    geom_errorbar(mapping = aes(x = name, ymin = value-se_dat_long$value, ymax = value+se_dat_long$value),
                  position = position_dodge(width = 0.9, preserve = "single"),
                  width = 0.5,
                  inherit.aes = TRUE) +
    scale_fill_viridis(discrete = TRUE, option = "G", begin = 0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8), aspect.ratio = 0.6) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, upper_lim))
  #theme_minimal() +
  #theme(axis.text.x = element_text(angle = 25, hjust = 0.8, size = 8))
}

xlab("Cell Type") + ylab("dMFI") + theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8), aspect.ratio = 0.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8500)) +
  scale_shape_manual(values = c(16,17,15,18))

bar_plot_fun("Analysis Results/cd8_pha_MFI_table.csv", factors, grps_of_interest)

bar_plot_fun("Analysis Results/cd8_pha_MFI_table.csv", factors, sex_four_group,upper_lim = 3)

bar_plot_fun("Analysis Results/cd8_pha_dMFI_linear_table.csv", factors, sex_four_group) +
  xlab("Cell Type") + ylab("dMFI") + theme(aspect.ratio = 0.6)



# ggplot(long, aes(x = name, y = value, fill = Group)) +
#   stat_summary(fun = median, geom = "bar", position = "dodge", color = "black") +
#   stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.5) +
#   scale_fill_viridis(discrete = TRUE, option = "G", begin = 0.2) +
#   geom_point(mapping = aes(name, value, shape = Group, fill = Group),
#              position = position_jitterdodge(jitter.width = 0.05, dodge.width = 1, seed = 42),
#              size = 1) +
#   theme_minimal()
