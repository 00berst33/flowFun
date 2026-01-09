library(openCyto)
library(CytoML)
library(flowGate)
library(data.table)

############
# openCyto #
############


# Set the directory you would like to store analysis results in.
work_dir <- file.path("C:/Users/00ber/OneDrive/Desktop/VPC") ### EDIT

## Load in raw data for first time
# The name of the directory containing the .fcs files to analyze
data_dir <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw" ### EDIT

# Get all filenames and create GatingSet
files <- list.files(data_dir, full.names = TRUE)
cs <- flowWorkspace::load_cytoset_from_fcs(files)
gs <- flowWorkspace::GatingSet(cs)

# Note: to clear, flowWorkspace::cs_cleanup_temp() and flowWorkspace::gs_cleanup_temp()

# Save resulting GatingSet to disk so that we have a backup
dir.create(file.path(work_dir, "backup_gs"))
flowWorkspace::save_gs(gs, path = file.path(work_dir, "backup_gs"))

# NOTE: if you would like to discard the edits you have made to the current
# GatingSet object, reload the backup from disk an make another deep copy, then
# proceed.
gs1 <- flowWorkspace::load_gs(file.path(work_dir, "backup_gs"))
gs1 <- flowWorkspace::gs_clone(gs1)

# Specify L/D stain
ld_stain <- "BUV496-A"

# Gating template csv
# Try editing the data.table yourself if you would like to add more gates.
# For explanations on how to do this, try running each command in the R console:
#   `vignette(package="openCyto", "openCytoVignette")`
#   `vignette(package="openCyto", "HowToAutoGating")`
# They each will take you to page that shows you how this file should look, and
# what it can be used for.
gt_table <- data.table(alias = c("nonMargins", "nonDebris", "singlets", "live"),
                       pop = c("+", "+", "+", "-"),
                       parent = c("root", "nonMargins", "nonDebris", "singlets"),
                       dims = c("FSC-A,SSC-A", "FSC-A", "FSC-A,FSC-H", "BUV496-A"),
                       gating_method = c("boundary", "gate_mindensity", "singletGate", "gate_mindensity"),
                       gating_args = c("min=c(0,0),max=c(262143,262143)", "gate_range=c(15000,80000)", NA, "gate_range=c(1,2)"),
                       collapseDataForGating = c(TRUE, TRUE, TRUE, TRUE),
                       groupBy = c(46, 46, 46, 46),
                       preprocessing_method = c(NA, NA, NA, NA),
                       preprocessing_args = c(NA, NA, NA, NA))
gt <- openCyto::gatingTemplate(gt_table)

# Compensation matrix csv
comp_mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/morgans_comp_matrix.csv", ### EDIT: path to compensation matrix
                     check.names = FALSE) # note check.names


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
gt_gating(gt, gs1)

###
# This is where I have begun to get the following error: "Invalid input type, expected 'double' actual 'list'".
# If this happens try commenting out the lines with `collapseDataForGating` and `groupBy` and running everything after again.
# Also try running lines 141 and 153 with transformed and compensated GatingSet

# ^ fixed by deleting all related data and redo-ing analysis


## Check one gate for each sample
# nonDebris gate
plotAllSamples(gs1, "FSC-A", "SSC-A", "nonMargins", "nonDebris")

# singlets gate
plotAllSamples(gs1, "FSC-A", "FSC-H", "nonDebris", "singlets")

# live cell gate
plotAllSamples(gs1, "BUV496-A", "FSC-A", "singlets", "live")

##### Redraw gates manually, if desired
### For one sample
# Check gate before
autoplot(gs1[[7]])

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
               dims = c("BUV496-A", "FSC-A"),
               ref_sample = 42) # sample to draw gate on; applied to all other samples when no `sample_ids` given

# live cell gate after
plotAllSamples(gs1, "BUV496-A", "FSC-A", "singlets", "live")

# save and/or export here


#####
# Add sequentially instead of loading default CSV
# May need to make fresh clone beforehand

# number of samples
num_samples <- length(gs1)

# nonMargins
gs_add_gating_method(gs1,
                     alias = "nonMargins",
                     pop = "+",
                     parent = "root",
                     dims = "FSC-A,SSC-A",
                     gating_method = "boundary",
                     gating_args = "min=c(0,0),max=c(262143,262143)",
                     collapseDataForGating = TRUE,
                     groupBy = 1
                     ) # find one gate for all samples


# nonDebris
gs_add_gating_method(gs1,
                     alias = "nonDebris",
                     pop = "+",
                     parent = "nonMargins",
                     dims = "FSC-A",
                     gating_method = "gate_mindensity",
                     #gating_args = "gate_range=c(10000,80000)",
                     collapseDataForGating = TRUE,
                     groupBy = 1
                     )

# Plot debris gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `SSC-A`), subset = "nonMargins") + #add subset attribute before here?
  geom_hex(bins = 50) +
  geom_gate("nonMargins") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)


# singlets
gs_add_gating_method(gs1,
                     alias = "singlets",
                     pop = "+",
                     parent = "nonDebris",
                     dims = "FSC-A,FSC-H",
                     gating_method = "singletGate",
                     collapseDataForGating = TRUE,
                     groupBy = num_samples)

# Plot singlet gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `FSC-H`), subset = "nonDebris") +
  geom_hex(bins = 50) +
  geom_gate("singlets") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)


# live cells
gs_add_gating_method(gs1,
                     alias = "live",
                     pop = "-",
                     parent = "singlets",
                     dims = "BUV496-A",
                     gating_method = "gate_mindensity",
                     gating_args = "gate_range=c(1,2)",
                     collapseDataForGating = TRUE,
                     groupBy = num_samples)

# Live cell gate
ggcyto::ggcyto(gs1, mapping = aes(x = !!enquo(ld_stain), y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 50) +
  geom_gate("live") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)


gt_df <- data.frame(alias = c(),
                    pop = c(),
                    parent = c(),
                    dims = c(),
                    gating_method = c(),
                    gating_args = c(),
                    collapseDataForGating = c(),
                    groupBy = c(),
                    preprocessing_method = c(),
                    preprocessing_args = c())


###
# Plots for all samples

# Change visibility of root node
# gs_pop_set_visibility(gs1, "root", FALSE)

# Plot debris gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `SSC-A`), subset = "nonMargins") + #add subset attribute before here?
  geom_hex(bins = 70) +
  geom_gate("nonDebris") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Singlet gate
ggcyto::ggcyto(gs1, mapping = aes(x = `FSC-A`, y = `FSC-H`), subset = "nonDebris") +
  geom_hex(bins = 70) +
  geom_gate("singlets") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)

# Live cell gate
ggcyto::ggcyto(gs1, mapping = aes(x = !!enquo(ld_stain), y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 70) +
  geom_gate("live") +
  theme(text = element_text(size = 4)) +
  geom_stats(size = 1)



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


# Must download docker image, see ?gatingset_to_flowjo()
# To save current GatingSet as a FlowJo workspace;
CytoML::gatingset_to_flowjo(gs1, file.path(getwd(), "test_human_gs.wsp"))
