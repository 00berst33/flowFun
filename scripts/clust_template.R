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
  # notes on how to determine for users
cols_to_cluster <- c(12, 14:16, 18, 20:25, 27:34, 36)

## Perform clustering
# Returns FlowSOM object
  # cytoset/GatingSet to flowSet
fsom <- FlowSOM::FlowSOM(fs,
                         colsToUse = cols_to_cluster,
                         seed = 42)

# Returns tidytable
  # remove?
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
