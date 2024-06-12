# Perform preprocessing of flow cytometry setup files. The steps taken are as
# follows: 
#   1. Removing margin events (i.e. values that fall outside the range the 
#      flow cytometer should be able to pick up).
#   2. Removing doublets.
#   3. Removing debris.
#   4. Performing a log-icle transformation on each channel.
# 
# A .pdf file plotting the results of removing doublets and debris from each file
# will be returned to check for any mistakes. There is also an option to apply
# a compensation matrix to your setup files and create a .pdf where, for each
# file, the channel/marker it corresponds to is plotted against every other
# channel/marker, so that you may confirm that the compensation matrix is appropriate.
# 
# This script requires the packages ggplot2, gridExtra, ggpointdensity, and viridis
# for plotting, and flowCore, PeacoQC, and CytoExploreR for preprocessing.


# Load the libraries necessary to run this script.
library(ggplot2)
library(gridExtra)
library(flowCore)
#library(FlowSOM)
library(PeacoQC)
library(CytoExploreR)

# Note that this workflow creates a directory structure as follows:
# 
# home ---------- Data ----- SetUp
#           +            +-- Preprocessed SetUp
#           +            +-- Raw
#           +            +-- Preprocessed Raw
#           +            +-- Aggregates
#           +            +-- Clustered Raw
#           +            +-- directories for any specific cell types of interest (e.g. CD4 T cells)
#           +
#           + --- Info
#           +
#           + --- RDS ------ Unedited
#           +            +-- Edited
#           +
#           + --- Preprocessing Results ----- QC plots
#           +
#           + --- Analysis Results ----- limma Analysis ----- Tables
#                                    +-- edgeR Analysis ----- Tables

# In this file, we are only concerned with the directories defined below. You 
# may change their names if you prefer, but do not change the names of the 
# variables. If the directories have already been created, and you would like 
# to rename them, you must do so in your file editor, then edit the 
# variables below accordingly.

# Will contain all .fcs files, including setup, raw, clustered, etc.
dir_data = "Data/"

# Contains raw setup files.
dir_setup = "Data/SetUp/"

# Contains preprocessed setup files.
dir_setup_prepr = "Data/Preprocessed SetUp/"

# Contains plots showing the results of our preprocessing steps, like removing 
# doublets and debris.
dir_prepr_res = "Preprocessing Results/"

# Checks whether or not the above directories exist, and creates them if not.
for (dir in c(dir_data, dir_setup, dir_prepr_res, dir_setup_prepr)) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Define a path to a reference setup file. The choice of file is arbitrary,
# we will use it to define the appropriate gates and transformations that will
# then be applied to all other setup files. 
ff_setup = read.FCS("Data/SetUp/SetUp_APC CD16.fcs")

# Define vector for renaming markers in setup files.
descs = c(NA, NA, NA, NA, NA, NA, "PHA-L", NA, NA, "CD123", NA, "CD16", "CD8", "CD3", 
          "IL10R", "HLA-DR", NA, "CD14", "CD11c", "TIM3", "CCR7", "CD38", "CD45RA", 
          "LD", "TCRgd", "CD56", "CD19", "CD39", "CD4", "PD1", "IgD", "CD25", NA, 
          "Fas", NA)

# Rename markers in file.
pData(parameters(ff_setup))$desc = descs

# Get all channel names, excluding time, FSC, and SSC.
all_channels = as.vector(pData(parameters(ff_setup))$name)[7:(nrow(parameters(ff_setup))-1)]

# Define helper functions for getting a marker's corresponding channel, or vice versa.
get_marker = function(ff, channel) {
  nm = as(pData(parameters(ff))[which(pData(parameters(ff))[, "name"] == channel),
                                "desc"], "character")
  if (is.na(nm)) {
    return(channel)
  } else {
    return(nm)
  }
}

get_channel = function(ff, marker) {
  return(as(pData(parameters(ff))[which(pData(parameters(ff))[, "desc"] == marker),
                                  "name"], "character"))
}

# Define a function for plotting events that were removed between steps. 
plot_before_after = function(ff1, ff2, channel1, channel2, ncells) {
  df = data.frame(x = exprs(ff1)[, channel1],
                  y = exprs(ff1)[, channel2])
  i = sample(nrow(df), ncells)
  if (!"Original_ID" %in% colnames(exprs(ff1))) {
    ff1@exprs = cbind(ff1@exprs,
                      Original_ID = seq_len(nrow(ff1@exprs)))
  }
  plot = ggplot(df[i,], aes(x = x, y = y)) + 
    geom_point(size = 0.5,
               color = ifelse(exprs(ff1)[i,"Original_ID"] %in% 
                                exprs(ff2)[,"Original_ID"], 'blue', 'red')) + 
    xlab(get_marker(ff1, channel1)) + 
    ylab(get_marker(ff1, channel2)) +
    theme_minimal() + theme(legend.position = "none")
  ggpubr::ggarrange(plot)
}

# Remove margin events.
ff_setup_m = RemoveMargins(ff_setup, c("FSC-A", all_channels))

# Remove doublets.
ff_setup_d = RemoveDoublets(ff_setup_m)

# Check results of removing doublets. Removed cells are hightlighted in red.
plot_before_after(ff_setup_m, ff_setup_d, "FSC-A", "FSC-H", 10000)

# Draw a gate to remove debris. A window will pop up on your computer;
# you may draw a polygon gate normally as you would in FlowJo, but to close it,
# instead of double-clicking, click "Stop", then "Stop locator" in the top left 
# corner, and close the window.
ff_setup_gate = cyto_gate_draw(x = ff_setup_d,
                               alias = "non debris",
                               channels = c("FSC-A", "SSC-A"))
ff_setup_gate = ff_setup_gate$`non debris`

# Apply previously defined gate to data.
ff_setup_g = ff_setup_d[filter(ff_setup_d, ff_setup_gate)@subSet, ]

# Check results of removing debris from the data.
plot_before_after(ff_setup_d, ff_setup_g, "FSC-A", "SSC-A", 10000)

# Define a logicle transformation on all channels except time, FSC, and SSC.
trans_setup = estimateLogicle(ff_setup_g, all_channels)

# Apply the transformation to all channels.
ff_setup_t = transform(ff_setup_g, trans_setup)

# Get all setup files that require preprocessing.
setup_files = list.files(path = dir_setup,
                         pattern = "*.fcs")

# The following loop applies the gates and transformation defined above to all
# setup files, and returns the results in a .pdf file hightlighting the removed
# events for each sample.
pdf(paste0(dir_prepr_res, "setup_preprocessing.pdf"))
for (file in setup_files) {
  ff_setup = read.FCS(paste0(dir_setup, file))
  pData(parameters(ff_setup))$desc = descs
  ff_setup_m = RemoveMargins(ff_setup, c("FSC-A", all_channels))
  ff_setup_d = RemoveDoublets(ff_setup_m)
  
  p1 = plot_before_after(ff_setup_m, ff_setup_d, "FSC-A", "FSC-H", 5000)
  
  ff_setup_g = ff_setup_d[filter(ff_setup_d, ff_setup_gate)@subSet, ]
  
  p2 = plot_before_after(ff_setup_d, ff_setup_g, "FSC-A", "SSC-A", 5000)
  
  grid.arrange(p1, p2, ncol = 2, nrow = 1, top = file)
  
  write.FCS(ff_setup_g,
            file = paste0(dir_setup_prepr, file))
  
  ff_setup_t = transform(ff_setup_g, trans_setup)
  
  write.FCS(ff_setup_t,
            file = paste0(dir_setup_prepr_trans, file))
}
dev.off()


###############
# The following code plots the results of compensating our setup data.
library(ggpointdensity)
library(viridis)

# Define a new function for creating each plot in the .pdf file.
ggplot_cells = function(ff, channel1, channel2) {
  df = data.frame(x = exprs(ff)[, channel1],
                  y = exprs(ff)[, channel2])
  i = sample(nrow(df), 5000)
  plot = ggplot(df[i,], aes(x = x, y = y)) + 
    geom_pointdensity(size = 0.01, show.legend = FALSE) +
    scale_color_viridis(option = "turbo") +
    xlab(GetMarkers(ff, channel1)) + 
    ylab(GetMarkers(ff, channel2)) + 
    theme(panel.background = element_blank())
  return(plot)
}

# Get names of all preprocessed setup files.
files = list.files(path = dir_setup_prepr,
                   pattern = ".fcs",
                   full.names = FALSE)

# Get all channel names, excluding time, FSC, and SSC.
all_channels = as.vector(pData(parameters(ff_setup))$name)[7:(nrow(parameters(ff_setup))-1)]

#setup_channels = as.vector(GetChannels(ff_setup, as.vector(na.exclude((pData(parameters(ff_setup))$desc)))))

# Read in compensation matrix.
# Assumes matrix is in working directory.
comp_matrix = read.csv("morgans_comp_matrix.csv", check.names = FALSE)

auto_comp_matrix = read.csv("auto_final.csv", check.names = FALSE)[-1]
colnames(auto_comp_matrix) = sub(" ::.*", "", colnames(auto_comp_matrix))

# Create .pdf file to check results of compensation. For each file, the 
# channel/marker it corresponds to is plotted against every other channel.
# Note that the loop takes a couple minutes to run.
pdf("compensation_results.pdf")
for (file in files) {
  print(paste("Processing", file))
  ff = read.FCS(file.path("Data/Preprocessed SetUp", file))
  
  if (file == "SetUp_Unstain.fcs") {
    file_channel = "FITC-A"
  } else if (file == "SetUp_BV421 CD45RA.fcs") {
    file_channel = "BV421-A"
  } else {
    file_channel = GetChannels(ff, sub("^SetUp_(?:[^ ]+ )?([^\\.]+)\\.fcs$", "\\1", file))
  }
  
  ff_c = compensate(ff, auto_comp_matrix)
  
  trans_new = estimateLogicle(ff_c, all_channels)
  ff_t = transform(ff_c, trans_new)
  
  plots = lapply(1:length(setup_channels), function(i)
    ggplot_cells(ff_t, file_channel, setup_channels[i]))
  
  grid.arrange(grobs = plots, top = file, newpage = TRUE)
}
dev.off()

