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
#library(CytoML)
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

###
## OR, Load in existing data
# note: changes made to object created from load_cytoset() do not affect
#   original data unless `backend_readonly=FALSE`
# Load previously created GatingSet, if needed
gs <- flowWorkspace::load_gs(file.path(work_dir, "backup_gs"))
# cs <- flowWorkspace::load_cytoset(file.path(work_dir, "raw_cs"))
###

# Make deep clone of GatingSet
# flowWorkspace::gs_cleanup_temp(gs1) # run first if overwriting old `gs1` edited with doPreprocessing()
gs1 <- flowWorkspace::gs_clone(gs)

# Compensate
comp_mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/morgans_comp_matrix.csv",
                     check.names = FALSE) # note check.names

# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.
# Check
colnames(comp_mat)
colnames(gs)

# Standardize marker/channel names if necessary
# colnames(comp_mat) <-
# cs_swap_colnames(), cs_rename_channel(), cs_rename_marker()

# Define custom transformation, if desired
# Example:
# trans <- flowCore::estimateLogicle(gs[[1]],
#                                    channels = colnames(comp_mat))

# Specify stain for dead cells, if you have one
ld_stain <- "BUV496-A"

# Preview preprocessing steps
# !!! Note: debris_args, singlet_args, and live_args must be strings, specifying arguments
# that are passed to openCyto::gate_mindensity and/or flowStats::gate_singlet
previewPreprocessing.GatingSet(gs1,
                               sample_ind = 1,
                               compensation = comp_mat,
                               transformation = "logicle", # may be a custom `transformList`
                               ld_channel = ld_stain,
                               debris_args = "gate_range = c(0, 75000)",
                               singlet_args = list(),
                               live_args = list())

# Apply preprocessing steps to all samples
doPreprocessing.GatingSet(gs1,
                          compensation = comp_mat,
                          transformation = "logicle",
                          ld_channel = ld_stain,
                          debris_args = "gate_range = c(5000, 75000)",             # see also `openCyto::gate_mindensity()`
                          singlet_args = list(),
                          live_args = "gate_range = c(0, 2)")                   # see also `flowStats::gate_singlet()`

## Make graphs to check results of preprocessing:
# Non-debris gate

attributes(gs1)$subsets <- # add parent with margin events removed here?

ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `SSC-A`), subset = "nonMargins") + #add subset attribute before here?
  geom_hex(bins = 70) +
  geom_gate("nonDebris") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Singlet gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `FSC-H`), subset = "nonDebris") +
  geom_hex(bins = 100) +
  geom_gate("singlets") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Live cell gate
ggcyto::ggcyto(gs1, mapping = aes(x = !!enquo(ld_stain), y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 100) +
  geom_gate("live") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Once you have results, save and convert data to flowSet to perform clustering
flowWorkspace::save_gs(gs1, path = file.path(work_dir, "preprocessed_gs")) ### should cytoset be saved instead?
fs <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs1))

############
# openCyto #
############

# Load gs on disk
gs1 <- flowWorkspace::load_gs(file.path(work_dir, "backup_gs"))
gs1 <- flowWorkspace::gs_clone(gs1)

# if you are adding gates sequentially instead of all in one .csv,
# clear existing gates with gs_add_gating_method_init()

# Gating template csv
# Try editing the CSV yourself if you like, running:
#   `vignette(package="openCyto", "openCytoVignette")`
#   `vignette(package="openCyto", "HowToAutoGating")`
# will take you to page that shows you how this file should look and how to edit it appropriately.
gt_file <- "C:/Users/00ber/Downloads/gtfile_basic.csv"
gt <- openCyto::gatingTemplate("C:/Users/00ber/Downloads/gtfile_basic.csv")

# Compensation matrix csv
comp_mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/morgans_comp_matrix.csv",
                     check.names = FALSE) # note check.names


# Visualize current gating scheme
plot(gt)

# custom function
gateMargins <- function(ff, channels = colnames(flowCore::exprs(ff))) {
  ran <- Biobase::pData(flowCore::parameters(ff)) %>%
    dplyr::select(name, minRange, maxRange) %>%
    dplyr::filter(name %in% channels)
  rownames(ran) <- ran$name
  bounds <- ran %>% dplyr::select(-name)

  gate <- flowCore::rectangleGate(t(bounds), filterId = "nonMargins")
  return(gate)
}

dat = flowWorkspace::gs_pop_get_data(gs1)

for (sn in sampleNames(gs1)) {
  ff <- flowWorkspace::gh_pop_get_data(gs1[[sn]], "root")   # get flowFrame for root node
  ff <- flowWorkspace::cytoframe_to_flowFrame(ff)
  gate <- gateMargins(ff)
  flowWorkspace::gs_pop_add(gs1, gate, parent = "root")    # add gate under root population
  flowWorkspace::recompute(gs1, "nonMargins", sn)   # recompute population for this sample
}

# register
openCyto::register_plugins(gateMargins, methodName = "gateMargins", dep = "flowCore", "gating")


# Perform quality control with PeacoQC
for (i in seq_along(gs1)) {
  ff <- flowWorkspace::gh_pop_get_data(gs1[[i]])
  ff <- flowWorkspace::cytoframe_to_flowFrame(ff)
  qc_ff <- PeacoQC::RemoveMargins(ff,
                                  channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H"))

  # Replace data inside GatingSet
  qc_ff <- flowWorkspace::flowFrame_to_cytoframe(qc_ff)
  identifier(qc_ff) <- identifier(ff)
  gs1[[i]] <- qc_ff
}


# Compensate
# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.
# Check
colnames(comp_mat)
colnames(gs1)

compensate(gs1, comp_mat)


# Transform
trans <- flowCore::estimateLogicle(gs1[[1]], channels = colnames(comp_mat))
flowCore::transform(gs1, trans)

pp_gs1 <- flowWorkspace::gs_clone(gs1)


# Apply gating scheme
gt_gating(gt, gs1)


# Check gates
autoplot(gs1[[1]], bins = 100, axis_inverse_trans = FALSE)

# Redraw gates if desired
# gh <- gs1[[ind]]
# flowWorkspace::gh_pop_set_gate(gh, "/nonMargins/nonDebris", new_filter)
# flowWorkspace::recompute(gh, "/nonMargins/nonDebris")

# Save final gating
flowWorkspace::save_gs(gs1, "preprocessed_gs")


# preprocessing controls
#####
## Load in raw data for first time
# The name of the directory containing the .fcs files to analyze
fmo_dir <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/FMO" ###

# Get all filenames and create GatingSet
fmo_files <- list.files(fmo_dir, full.names = TRUE)
fmo_cs <- flowWorkspace::load_cytoset_from_fcs(fmo_files)
fmo_gs <- flowWorkspace::GatingSet(fmo_cs)

compensate(fmo_gs, comp_mat)

# Transform
# trans <- flowCore::estimateLogicle(gs1[[1]], channels = colnames(comp_mat))
flowCore::transform(fmo_gs, trans)


# Apply gating scheme
gt_gating(gt, fmo_gs)

# Save
save_gs(fmo_gs, "fmo_gs")

# Or, if some gates were manually readjusted for raw data and you would like
# to carry those over:

# copyTree(gs1, gs_fmo)



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


# After clustering; how to save results? write to cytoset? how to gate on subsets?

##############
# Clustering #
##############

# clustering GatingSet
ex_fs <- flowWorkspace::gs_pop_get_data(gs1, "live")
ex_fs <- flowWorkspace::cytoset_to_flowSet(ex_fs)
agg <- FlowSOM::AggregateFlowFrames(ex_fs, cTotal=50000, writeOutput=FALSE)
fsom <- FlowSOM::FlowSOM(agg,
                         xdim = 10,
                         ydim = 10,
                         nClus = 10)
saveFlowSOMEmbedding(fsom, "fsom_minimal.rds")

# for controls
# clustering GatingSet
ex_fs <- flowWorkspace::gs_pop_get_data(fmo_gs, "live")
ex_fs <- flowWorkspace::cytoset_to_flowSet(ex_fs)
fmo_agg <- FlowSOM::AggregateFlowFrames(ex_fs, cTotal=50000, writeOutput=FALSE)
fmo_clus <- loadFlowSOMEmbedding("fsom_minimal.rds", fmo_agg) # use saved FlowSOM embedding to cluster


#############
# DA and DE #
#############

# Read in existing clustering (or use clustering created above)
fsom_26 <- readRDS("C:/Users/00ber/OneDrive/Desktop/VPC/human1/RDS/Edited/fsom_edited_26_meta.rds")

# Save mapping of the clustering you are happy with
saveFlowSOMEmbedding(fsom_26, "fsom26_minimal.rds")

# Change GatingSet to flowSet
fs <- flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(gs1))

# Subset samples to a less computationally intense number of cells (`cTotal` may be edited)
fs <- FlowSOM::AggregateFlowFrames(fs, cTotal = 250000)

# Map this subset to the given clustering
new_fsom <- loadFlowSOMEmbedding("fsom26_minimal.rds", fs)

# Convert resulting FlowSOM object to a data.table
fsom_dt <- flowSOMToTable(new_fsom)

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
                                 filename_col = "X...File.Name",
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

# Change scale and apply controls here, if desired
# ...

# Set markers of interest
marker_cols <- c("BV711-A", "FITC-A")

# Find expression matrix: metacluster.marker by sample
collapsed <- getExprMatDE(fsom_dt, marker_cols)

# Create linear models.
# NOTE: weights questionable
lm_model <- limma::lmFit(object = collapsed,
                         design = design)

# Perform statistical tests.
contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)

# View results
tt = limma::topTable(limma_ebayes)



# ex_k <- fsApply(gs_pop_get_data(gs1,y="live"),flowCore::exprs)
# agg <- FlowSOM::AggregateFlowFrames(ex_k, cTotal=50000, writeOutput=FALSE)
# can this be added to gs as boolean gate?
# flowWorkspace::flowFrame_to_cytoframe() # then save

# ex_mats <- flowWorkspace::gs_get_singlecell_expression_by_gate(
#   gs1,
#   nodes = "live",
#   swap = TRUE,
#   other.markers = colnames(comp_mat),
#   threshold = FALSE
#   ) # use if only interested in doing DE and DA for previously gated data?

###
# TO DO #
###
  # add column of metacluster to FMO control data matrices
    # only original clustering and each related subclustering needs to be saved

#####

# ignore probably
# # normally you are finding matrix where rows are samples and columns are metaclusters:
fsom_dt <- fsom_dt %>%
  dplyr::select(c(`FITC-A`, File, Metacluster)) %>%
  dplyr::group_by(File, Metacluster) %>%
  dplyr::summarise(med = median(`FITC-A`)) %>%
  tidyr::pivot_wider(names_from = File, values_from = med)

rownames(fsom_dt) <- fsom_dt$Metacluster
fsom_dt <- fsom_dt %>%
  dplyr::select(-Metacluster)
# do this for each marker we wish to pass to limma
# you are making lists where each item is a marker, but instead each item may be a metacluster instead,
# thus rows are markers instead of metaclusters for each item

# if you want it back-transformed, do that before passing it to function


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
