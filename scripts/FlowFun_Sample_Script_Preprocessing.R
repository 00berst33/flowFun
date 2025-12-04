# Load packages.
library(flowCore)
library(flowFun)
library(flowWorkspace)
library(openCyto)
library(CytoML)
library(flowGate)
library(data.table)
library(ggcyto)
library(ggplot2)
library(flowDensity) #from workflowVignetteUpdated

############
# openCyto #
############
## This part of the code takes you from raw fcs files to compensated data with only live cells included (or whatever first gating sets you are interested in)

# Set the directory you would like to store analysis results in.
work_dir <- file.path("...") ##ADD IN YOUR PATHNAME

#####
## Load in raw data for first time
# The name of the directory containing the .fcs files to analyze
data_dir <- "..." ##ADD IN YOUR PATHNAME

# Get all filenames and create GatingSet
files <- list.files(data_dir, full.names = TRUE)
cs <- flowWorkspace::load_cytoset_from_fcs(files)
gs <- flowWorkspace::GatingSet(cs)

# Note: to clear, flowWorkspace::cs_cleanup_temp() and flowWorkspace::gs_cleanup_temp() ### MR: What is this?

# Save resulting GatingSet to disk so that we have a backup ### MR: Did we decide to remove this? I took this from
flowWorkspace::save_gs(gs, path = file.path(work_dir, "backup_gs"))

# NOTE: if you would like to discard the edits you have made to the current
# GatingSet object, reload the backup from disk an make another deep copy, then
# proceed.
gs1 <- flowWorkspace::load_gs(file.path(work_dir, "backup_gs"))
gs1 <- flowWorkspace::gs_clone(gs1)

# Specify viability stain
ld_stain <- "BUV496-A"

# Gating template csv
# Try editing the data.table yourself if you would like to add more gates.
# For explanations on how to do this, try running each command in the R console:
#   `vignette(package="openCyto", "openCytoVignette")`
#   `vignette(package="openCyto", "HowToAutoGating")`
# They each will take you to page that shows you how this file should look, and
# what it can be used for.

num_samples <- length(gs1)

gt_table <- data.table(alias = c("nonMargins", "nonDebris", "singlets", "live"),
                       pop = c("+", "+", "+", "-"),
                       parent = c("root", "nonMargins", "nonDebris", "singlets"),
                       dims = c("FSC-A,SSC-A", "FSC-A", "FSC-A,FSC-H", ld_stain),
                       gating_method = c("boundary", "gate_mindensity", "singletGate", "gate_mindensity"),
                       gating_args = c("min=c(0,0),max=c(262143,262143)", "gate_range=c(15000,80000)", NA, "gate_range=c(1,2)"),
                       collapseDataForGating = c(TRUE, TRUE, TRUE, TRUE),
                       groupBy = seq(num_samples, 4),
                       preprocessing_method = c(NA, NA, NA, NA),
                       preprocessing_args = c(NA, NA, NA, NA))
gt <- openCyto::gatingTemplate(gt_table)

# Compensation matrix csv
comp_mat <- read.csv("...", ### EDIT: path to compensation matrix
                     check.names = FALSE) # note check.names ### MR: might be good to say why


# Visualize current gating scheme
plot(gt)


# Compensate
# !!! Note: The column names of your compensation matrix must match those
#     found in your GatingSet.
#   Check that that is the case here:
colnames(comp_mat)
colnames(gs1)

compensate(gs1, comp_mat)

# Transform
trans <- flowCore::estimateLogicle(gs1[[1]], channels = colnames(comp_mat))
flowCore::transform(gs1, trans)


# Apply gating scheme
### NOTE: Make sure you do not have tidytable loaded before you run this, if you do, you might get the following error: "Invalid input type, expected 'double' actual 'list'".
gt_gating(gt, gs1)

## Check each gate across all the samples  ### MR: Note that if you have many samples it's hard to see, would be better to include info for checking a subset
# nonDebris gate
plotAllSamples(gs1, "FSC-A", "SSC-A", "nonMargins", "nonDebris")

# singlets gate
plotAllSamples(gs1, "FSC-A", "FSC-H", "nonDebris", "singlets")

# live cell gate
plotAllSamples(gs1, !!enquo(ld_stain), "FSC-A", "singlets", "live")

# IF YOU ARE HAPPY WITH THE GATES, run this line of code to save and then move on to clustering
# Otherwise, skip this line and proceed to next part for redrawing
flowWorkspace::save_gs(gs1, path = file.path(work_dir, "processed_data_gs"))


##### Redraw gates manually, if desired
### For one sample
# Check gate before
autoplot(gs1[[7]]) # 7 here denotes is the sample, you can plot different samples by changing the number
### MR: The bining is weird here (get big hexagons instead of individual dots)

# Redraw gate
editGateManual(gs1,
               node = "nonDebris",
               dims = c("FSC-A", "SSC-A"),
               sample_ids = 7)

# Check gate after
autoplot(gs1[[7]])

### For a subset of samples
# Redraw gate
editGateManual(gs1,
               node = "singlets",
               # dims = c("FSC-A", "FSC-H"), # if `dims` isn't specified, the dimensions used in original gate are used
               ref_sample = 25,
               sample_ids = 23:46) # sample_ids may also be a vector of sample indices/names to apply the new gate to

# singlets gate after
plotAllSamples(gs1, "FSC-A", "FSC-H", "nonDebris", "singlets")
# may also plot a subset of the GatingSet, for example:
plotAllSamples(gs1[23:46], "FSC-A", "FSC-H", "nonDebris", "singlets")

### For all samples
editGateManual(gs1,
               node = "live",
               dims = c(!!enquo(ld_stain), "FSC-A"),
               ref_sample = 42) # sample to draw gate on; applied to all other samples when no `sample_ids` given

# live cell gate after
plotAllSamples(gs1, !!enquo(ld_stain), "FSC-A", "singlets", "live")

### ONCE HAPPY WITH GATES SAVE THE PROCESSED DATA
flowWorkspace::save_gs(gs1, path = file.path(work_dir, "processed_data_gs"))
