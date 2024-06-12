# Perform preprocessing of raw flow cytometry data. The steps taken, in order,
# are: 
#   1. Removing margin events (i.e. values that fall outside the range the 
#      flow cytometer should be able to pick up).
#   2. Removing doublets.
#   3. Removing debris.
#   4. Compensating the data.
#   5. Performing a log-icle transformation on each channel. 
#   6. Removing dead cells. 
#   7. Quality control, to check and correct for occurrences like clogs, speed 
#      changes, etc. and flag any files that may be too low-quality to include
#      in the analysis.
#
# The script will output a few files to check for any mistakes. The first is
# a pdf highlighting the cells that were removed in each .fcs file. The others
# are plots for any files that were flagged during quality control, where each 
# channel is plotted against time and removed events are marked.
# 
# Finally, this script requires that the packages ggplot2, gridExtra, flowCore,
# flowDensity, flowCut, PeacoQC, and CytoExploreR have been installed.

# Load necessary packages.
library(ggplot2)
library(gridExtra)
library(flowCore)
library(flowDensity)
library(PeacoQC)
library(CytoExploreR)
library(flowCut)

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
#                                    +-- edgeR Analysis

# The purpose of each directory will be explained further as they become relevant.
# For preprocessing raw data, we are only concerned with the directories
# defined below. You may change their names if you prefer, but do not change the 
# names of the variables. If the directories have already been created, and you 
# would like to rename them, you must do so in your file editor, then edit the 
# variables below accordingly for the script to work properly.

# Contains all .fcs files organized in multiple sub-directories. 
dir_data = "Data/"

# Contains all raw .fcs files.
dir_raw = "Data/Raw/"

# Will contain the cleaned .fcs files.
dir_prepr = "Data/Preprocessed Raw2/"

# Will contain plots showing results of preprocessing.
dir_prepr_res = "Preprocessing Results/"

# Will contain the plots resulting from quality control.
dir_prepr_res_qc = "Preprocessing Results/QC plots"

# Define the name of a reference .fcs file, and read it in as a flowFrame object.
# The choice of file is arbitrary, we will use it to define the appropriate gates 
# and transformations that will then be applied to all other data files. 
ref_filename = "Ab_PHA_Ctrl AWB2.fcs"

# Define filepath to compensation matrix, which should be a .csv file. 
comp_matrix_filepath = "auto_final.csv"

# Define channel corresponding to live/dead cells. This should be a string.
ld_channel = "BUV496-A"

# A string specifying which plots you would like to be saved during quality
# control; either "All" or "Flagged Only".
Plot = "Flagged Only"

#################
# PREPROCESSING #
#################

# Checks whether the above directories already exist, and creates them if not.
# Note that at this point, it is expected that directories containing the raw
# data and setup data have already been created. 
for (dir in c(dir_data, dir_prepr, dir_prepr_res, dir_prepr_res_qc)) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

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

# Read in defined reference data file.
ff = read.FCS(paste0(dir_raw, ref_filename), truncate_max_value = FALSE)

# Get all channel names, excluding time, FSC, and SSC.
all_channels = as.vector(pData(parameters(ff))$name)[7:(nrow(parameters(ff))-1)]

# Remove margin events.
ff_m = RemoveMargins(ff, c("FSC-A", all_channels))

# Create density plot to check for doublets.
plotDens(ff_m, c("FSC-A", "FSC-H"))

# Remove doublets.
ff_d = RemoveDoublets(ff_m)

# Create a plot highlighting removed cells to check the result. The final
# parameter may be edited to include more or less cells in the plot, but note
# that runtime increases with the number of cells used.
plot_before_after(ff_m, ff_d, "FSC-A", "FSC-H", 10000)

# Similarly, create a density plot containing gated-on cells to check the results. 
plotDens(ff_d, c("FSC-A", "FSC-H"))

# Create a density plot to check for debris.
plotDens(ff_d, c("FSC-A", "SSC-A"))

# Ask user to draw a gate to remove debris via a pop-up window. NOTE: Drawing the 
# gate works similarly to FlowJo, where you click to place each point. However, to 
# close the gate, rather than double clicking, place your final point, and click 
# "Stop", then "Stop locator", in the top left of the window. The gate will then
# be saved to the "polygon_gate" variable and you may exit the pop-up window.
polygon_gate = cyto_gate_draw(x = ff_d,
                              alias = "non debris",
                              channels = c("FSC-A", "SSC-A"))
polygon_gate = polygon_gate$`non debris`

# Apply the gate to the data.
ff_g = ff_d[flowCore::filter(ff_d, polygon_gate)@subSet, ]

# Create a plot highlighting removed cells to check the result.
plot_before_after(ff_d, ff_g, "FSC-A", "SSC-A", 10000)

# Similarly, create a density plot to check the result. 
plotDens(ff_g, c("FSC-A", "SSC-A"))

# Read in compensation matrix, rename columns, and apply it to the data.
comp_matrix = read.csv(comp_matrix_filepath, check.names = FALSE)[-1]
colnames(comp_matrix) = sub(" :: .*", "", colnames(comp_matrix))
ff_c = compensate(ff_g, comp_matrix)

# Estimate a log-icle transformation on each channel and apply it to the data.
trans = estimateLogicle(ff_c, all_channels)
ff_t = transform(ff_c, trans)

# Create a plot to check for live/dead cells.
plotDens(ff_t, c(ld_channel, "FSC-A"))

# Ask the user to draw a gate on live cells.
ff_live_gate = cyto_gate_draw(x = ff_t,
                              alias = "live",
                              channels = c(ld_channel, "FSC-A"))
ff_live_gate = ff_live_gate$live

# Apply the live/dead gate to the data.
ff_l = ff_t[flowCore::filter(ff_t, ff_live_gate)@subSet, ]

# Create a plot highlighting removed cells to check the result.
plot_before_after(ff_t, ff_l, ld_channel, "FSC-A", 10000)

# Create a density plot to check results.
plotDens(ff_l, c(ld_channel, "FSC-A"))

#PQC = PeacoQC(ff = ff_l, 
#              channels = channels_of_interest,
#              determine_good_cells = "all", # may change to only MAD or IT
#              plot = TRUE,
#              save_fcs = FALSE,
#              output_directory = "Preprocessing Results/",
#              MAD = 6, # increase to make less strict
#              IT_limit = 0.55) # increase to make less strict

# Perform quality control using flowCut. Any files that:
#   1. are not monotonically increasing in time,
#   2. have sudden changes in fluorescence,
#   3. display a large gradual change of fluorescence in all channels,
#   4. display a very large gradual change of fluorescence in one channel,
# are flagged and have their relevant plots written to the 
# "Preprocessing Results/QC plots" directory. Edit the "Plot" variable
# above if you would like to change which plots are saved.
fc = flowCut(f = ff_l,
             FileID = ref_filename,
             Plot = Plot,
             Directory = dir_prepr_res_qc,
             Verbose = TRUE)
ff_fc = fc$frame

# Get the names of all raw .fcs files.
raw_files = list.files(path = dir_raw, 
                       pattern = ".fcs")

# Preprocess all raw data .fcs files, and write the cleaned files to the 
# "Data/Preprocessed Raw" directory. In addition, any generated plots are
# written within the "Preprocessing Results/" directory. It is likely that
# after the loop has finished running, multiple warnings will be given. These
# should be checked to see if any files are too low-quality to include in the
# analysis.
pdf(paste0(dir_prepr_res, "raw_preprocessing.pdf"))
for (file in raw_files) {
  print(file)
  ff = read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  ff_m = RemoveMargins(ff, c("FSC-A", all_channels))
  ff_d = RemoveDoublets(ff_m)
  
  p1 = plot_before_after(ff_m, ff_d, "FSC-A", "FSC-H", 5000)
  
  ff_g = ff_d[flowCore::filter(ff_d, polygon_gate)@subSet, ]
  
  p2 = plot_before_after(ff_d, ff_g, "FSC-A", "SSC-A", 5000)
  
  ff_c = compensate(ff_g, comp_matrix)
  ff_t = transform(ff_c, trans)
  ff_l = ff_t[flowCore::filter(ff_t, ff_live_gate)@subSet, ]
  
  p3 = plot_before_after(ff_t, ff_l, ld_channel, "FSC-A", 5000)
  
  grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = file)
  
  #PQC = PeacoQC(ff = ff_l, 
  #              channels = channels_of_interest,
  #              determine_good_cells = "all",
  #              plot = TRUE,
  #              save_fcs = FALSE,
  #              output_directory = "Preprocessing Results/",
  #              MAD = 6, 
  #              IT_limit = 0.55) 
  
  fc = flowCut(f = ff_l,
               FileID = file,
               Plot = Plot,
               Directory = dir_prepr_res_qc,
               Verbose = TRUE)
  ff_fc = fc$frame
  
  # Print a warning for any files with less than 60% live cells.
  if (nrow(ff_l)/nrow(ff_g) < 0.60) {
    warning(paste(file, "has less than 60% live cells."))
  }
  # Print a warning for any files that had more than 20% of its events removed by QC.
  if (nrow(ff_fc)/nrow(ff_l) < 0.80) {
    warning(paste(file, "had more than 20% of its events removed by flowCut."))
  }
  #if (PQC$PercentageRemoved > 20) { # CHANGE THRESHOLD AS DESIRED
  #  warning(paste(file, "had more than 20% of its events removed by PeacoQC quality control."))
  #}
  
  write.FCS(ff_fc,
            file = paste0(dir_prepr, file))
}
dev.off()