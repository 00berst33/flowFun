library(flowCore)
library(tidytable)
library(flowFun)

# Import the gating set (aka preprocessed FCS files) you generated with the preprocessing script
gs1 <- flowWorkspace::load_gs(file.path("..."))
gs1 <- flowWorkspace::gs_clone(gs1)

# clustering GatingSet
# Get cell population that made it through preprocessing
ex_fs <- flowWorkspace::gs_pop_get_data(gs1, "live")
# Convert data to flowSet so that it is compatible with FlowSOM
ex_fs <- flowWorkspace::cytoset_to_flowSet(ex_fs)
# Create aggregated .fcs file with a sample of cells to reduce computational load
agg <- FlowSOM::AggregateFlowFrames(ex_fs, cTotal=50000, writeOutput=FALSE)

# Check column numbers to figure out what to input for the next step (in colsToUse parameter)
data.table::as.data.table(list(cols = colnames(flowCore::exprs(agg))))

#OR if you want it more condensed use the following instead
# colnames(flowCore::exprs(agg))

fsom <- FlowSOM::FlowSOM(agg,
                         xdim = 12, #xdim and ydim together set the number of cluster (# = xdim x ydim, so here there will be 144)
                         ydim = 12,
                         colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34), #tells it which parameters to use for clustering
                         nClus = 25, # number of metclusters
                         transform = FALSE,
                         seed = 42) # make results reproducible
# Save clustering mapping so that it may be applied to similar datasets
saveFlowSOMEmbedding(fsom, "fsom_minimal.rds")

fsom_dt <- flowSOMToTable(fsom)

####################################################
# Exploring the Clusters and Reassigning as Needed #
####################################################

# In this part of the code you can explore the metaclusters and clusters to make sure you are happy
# happy with how everything is clustered. You don't have to follow the code top to bottom, instead
# use each visualization to probe the data and then can reassign clusters and/or merge metalclusters as needed.

### DEALING WITH METACLUSTERS ###

# GENERATE HEATMAP showing MFI for each clustering parameters for each metacluster
plotMetaclusterMFIs(fsom)

#if you want to just show specific parameters you can add c() with the parameters you want in the heatmap
# plotMetaclusterMFIs(fsom, c(14:16, 20))
# OR you can make a value that equals the parameters and input that
# cols_to_plot <- c(10, 12:14, 16, 18:23)
# plotMetaclusterMFIs(fsom, cols_to_plot)

# Generate 2D scatterplot for specific metaclusters
# In the plot, the clusters are indicated by the numbers and the metaclusters by the colouring
plotLabeled2DScatter(fsom_dt,
                     channelpair = colnames(flowCore::exprs(agg))[c(12,25)], #indicate the X and Y parameter you want to show
                     clusters = NULL,
                     metaclusters = 15) #indicate which metaclusters you want to plot
                                        #if you want to plot mutliple metacluses use c() instead if a number (e.g. c(1, 15))

# GENERATE UMAP coloured by metacluster.
### MR: The plot is cutting off the too close to the metacluster numbers
# num_cells is how many cells it includes in the UMAP
# seed is to make the UMAP consistent (same seed) or different (different seed) if you run it again
plotUMAP(fsom, num_cells = 5000, seed = 17)


#Merge metaclusters for FlowSOM object
fsom <- FlowSOM::UpdateMetaclusters(fsom,
                                    newLabels = c("9" = "1",
                                                  "2" = "1",
                                                  "8" = "4",
                                                  "23" = "15",
                                                  "21" = "15",
                                                  "20" = "15",
                                                  "22" = "15",
                                                  "7" = "6"))

# OR, if your data is a data.table made with flowSOMWrapper(), merge metaclusters with editTableMetaclusters()
# fsom_dt <- editTableMetaclusters(fsom_dt,
#                                  new_labels = c("9" = "1",
#                                                 "2" = "1",
#                                                 "8" = "4",
#                                                 "23" = "15",
#                                                 "21" = "15",
#                                                 "20" = "15",
#                                                 "22" = "15",
#                                                 "7" = "6"))


###### STOPPED HERE ######

# Merge metaclusters for datatable
#fsom <- FlowSOM::UpdateMetaclusters(fsom,
#                                 newLabels = c("9" = "1",
#                                               "2" = "1",
#                                              "8" = "4",
#                                             "23" = "15",
#                                            "21" = "15",
#                                           "20" = "15",
#                                          "22" = "15",
#                                         "7" = "6"))

# Label metaclusters ### MR: Need to fix for FlowSOM object
fsom <- FlowSOM::UpdateMetaclusters(fsom, newLabels = c("6" = "Monocytes",
                                                         "5" = "cDC",
                                                         "4" = "B cells",
                                                         "17" = "NK T cells",
                                                         "19" = "NK cells",
                                                         "1" = "CD8 T cells",
                                                         "11" = "gdT cells",
                                                         "15" = "CD4 T cells",
                                                         "18" = "pDC",
                                                         "3" = "Undefined"))


### DEALING WITH CLUSTERS ###

cols_to_cluster <- c(10, 12:14, 16, 18:23, 25:32, 34)

# Generate cluster MFI heatmap
plotClusterMFIs(fsom, cols_to_cluster, metaclusters = c(1, 25)) ### MR: The cluster numbers when you Zoom the plot are soooooo tiny...need to make bigger for old eyes like mine
#if you want to just show specific parameters you can replace cols_to_cluster with c() with the parameters you want in the heatmap
# plotClusterMFIs(fsom_dt, c(14:16, 20), metaclusters = c(3, 11))

# Generate 2D scatterplot
# In the plot, the clusters are indicated by the numbers and the metaclusters by the colouring
plotLabeled2DScatter(fsom_dt,
                     channelpair = c("APC-Cy7-A", "BUV563-A"),
                     clusters = c(4, 42, 76),
                     metaclusters = NULL)


# Merge metaclusters ### MR: What is this? Is it merging metaclusters or clusters?
# Both, `newLabels` merges metaclusters and `clusterAssignment` assigns clusters to a given metacluster
### MR: How could you reassign a cluster to a different existing metacluster instead of a new one (e.g. Undefined)
fsom <- FlowSOM::UpdateMetaclusters(fsom,
                                    newLabels = c("12" = "5",
                                                  "13" = "5",
                                                  "10" = "6",
                                                  "14" = "6",
                                                  "25" = "18",
                                                  "16" = "11",
                                                  "24" = "3"),
                                 clusterAssignment = c("4" = "Undefined",
                                                         "42" = "Undefined"))


##### ONCE YOU ARE HAPPY WITH YOUR CLUSTERS ####

# Label metaclusters
fsom <- FlowSOM::UpdateMetaclusters(fsom, newLabels = c("6" = "Monocytes",
                                                         "5" = "cDC",
                                                         "4" = "B cells",
                                                         "17" = "NK T cells",
                                                         "19" = "NK cells",
                                                         "1" = "CD8 T cells",
                                                         "11" = "gdT cells",
                                                         "15" = "CD4 T cells",
                                                         "18" = "pDC",
                                                         "3" = "Undefined"))

# Generate final heatmap
plotMetaclusterMFIs(fsom)

# Generate annotated heatmap
### MR: It would be good to have the original metacluster number on this plot
annotateMFIHeatmap(fsom)



################################
# Clustering a Subset of Cells #
################################

fsom_dt <- flowSOMToTable(fsom)

# If you want to define clusters of a subset of cells (i.e. one of the final metclusters from above),
# you can do that using the following code. An example of when to do this would be if you wanted to
# subset the CD8 T cells into various subtypes (e.g. naive, effector, memory, etc.)

# Get table of only CD8 T cells
new_table <- createFilteredAggregate(fsom_dt,
                                     num_cells = Inf,
                                     metaclusters = "CD8 T cells",
                                     clusters = NULL)

# In this example we are using principle components instead of defined parameters to do the clustering.
# This is helpful if you are not sure which markers to use for the clustering

# Perform PCA
pca_obj <- doPCA(new_table, cols_to_cluster)

# Draw scree plot
# This is to determine how many PCAs to include in the analysis. You are looking for an "elbow" in the plot
plotPCAScree(pca_obj)

# Recluster CD8 T cells
### How do we set the number of clusters for this?
fsom_sub <- clusterSubsetWithPCA(new_table,
                                 pca_obj = pca_obj,
                                 num_components = 6, # set this number based on where you see the elbow in the screen plot
                                 num_clus = 15, # this sets the number of clusters
                                 seed = 33)

#### YOU CAN NOW USE THE SAME FUNCTIONS AS ABOVE TO EXPLORE THE CLUSTERS ###
# some examples are included below


# Define new columns to use for the heatmaps - this is just the list of markers you want to see in the heatmaps
# you can choose them based on the cells you are looking at.
cols_of_interest <- c(15, 16, 18, 23:25, 28, 32)

# Generate heatmap for subsetted cells
plotMetaclusterMFIs(fsom_sub, cols_of_interest)
